"""Export Biomarker Intelligence Agent results to Markdown, JSON, PDF, CSV, and FHIR R4.

Provides public functions:
  - export_markdown()            -- human-readable report with evidence tables
  - export_json()                -- machine-readable structured data
  - export_pdf()                 -- styled PDF report via reportlab Platypus
  - export_csv()                 -- tabular CSV export for spreadsheet analysis
  - export_fhir_diagnostic_report() -- FHIR R4 DiagnosticReport JSON bundle

Author: Adam Jones
Date: March 2026
"""

import csv
import io
import json
import re
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from loguru import logger

from .models import (
    AnalysisResult,
    CrossCollectionResult,
    PatientProfile,
    SearchHit,
)

VERSION = "1.0.0"


# =====================================================================
# PUBLIC API
# =====================================================================


def generate_filename(extension: str) -> str:
    """Generate a timestamped filename with UUID suffix for export.

    Args:
        extension: File extension without dot (e.g. "md", "json")

    Returns:
        Filename like biomarker_report_20260301T143025Z_a1b2.md
    """
    import uuid
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = uuid.uuid4().hex[:4]
    return f"biomarker_report_{ts}_{suffix}.{extension}"


def export_markdown(
    query: str,
    response_text: str,
    evidence: Optional[CrossCollectionResult] = None,
    analysis: Optional[AnalysisResult] = None,
    filters_applied: Optional[dict] = None,
) -> str:
    """Export a query result as a Markdown report.

    Args:
        query: The user's original question.
        response_text: The LLM-generated response.
        evidence: CrossCollectionResult from retrieval.
        analysis: Optional AnalysisResult with patient analysis data.
        filters_applied: Dict of filters that were active.

    Returns:
        Complete Markdown report as a string.
    """
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    filters_str = _format_filters(filters_applied)

    lines = [
        "# Biomarker Intelligence Report",
        "",
        f"**Query:** {query}",
        f"**Generated:** {timestamp}",
        f"**Filters:** {filters_str}",
        "",
        "---",
        "",
        "## Response",
        "",
        response_text,
        "",
        "---",
        "",
    ]

    # Analysis summary section
    if analysis:
        lines.append("## Patient Analysis Summary")
        lines.append("")
        ba = analysis.biological_age
        lines.append(f"**Biological Age:** {ba.biological_age:.1f} years "
                      f"(chronological: {ba.chronological_age}, "
                      f"acceleration: {ba.age_acceleration:+.1f} years)")
        lines.append("")

        if analysis.critical_alerts:
            lines.append("### Critical Alerts")
            lines.append("")
            for alert in analysis.critical_alerts:
                lines.append(f"- {alert}")
            lines.append("")

        if analysis.disease_trajectories:
            lines.append("### Disease Risk Summary")
            lines.append("")
            lines.append("| Disease | Risk Level | Est. Years to Onset |")
            lines.append("|---------|-----------|-------------------|")
            for traj in analysis.disease_trajectories:
                years = f"~{traj.years_to_onset_estimate:.0f}" if traj.years_to_onset_estimate else "N/A"
                lines.append(f"| {traj.disease.value} | {traj.risk_level.value} | {years} |")
            lines.append("")

        if analysis.pgx_results:
            lines.append("### PGx Profile")
            lines.append("")
            for pgx in analysis.pgx_results:
                lines.append(f"- **{pgx.gene}** {pgx.star_alleles}: {pgx.phenotype.value}")
            lines.append("")

        lines.append("---")
        lines.append("")

    # Evidence section
    if evidence and evidence.hit_count > 0:
        lines.append("## Evidence Sources")
        lines.append("")
        lines.extend(_format_evidence_section(evidence))

        if evidence.knowledge_context:
            lines.append("")
            lines.append("## Knowledge Graph Context")
            lines.append("")
            lines.append(evidence.knowledge_context)

        # Search metrics
        lines.append("")
        lines.append("---")
        lines.append("")
        lines.append("## Search Metrics")
        lines.append("")
        lines.append("| Metric | Value |")
        lines.append("|--------|-------|")
        lines.append(f"| Total Results | {evidence.hit_count} |")
        lines.append(f"| Collections Searched | {evidence.total_collections_searched} |")
        lines.append(f"| Search Time | {evidence.search_time_ms:.0f}ms |")

    # Footer
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append(f"*Generated by HCLS AI Factory -- Biomarker Intelligence Agent v{VERSION}*")
    lines.append("")

    return "\n".join(lines)


def export_json(
    analysis_result: Optional[AnalysisResult] = None,
    query: str = "",
    response_text: str = "",
    evidence: Optional[CrossCollectionResult] = None,
) -> str:
    """Export analysis result as structured JSON.

    Args:
        analysis_result: AnalysisResult from patient analysis.
        query: Original question.
        response_text: LLM-generated response.
        evidence: CrossCollectionResult from retrieval.

    Returns:
        Pretty-printed JSON string.
    """
    timestamp = datetime.now(timezone.utc).isoformat()

    data: Dict[str, Any] = {
        "report_type": "precision_biomarker_analysis",
        "version": VERSION,
        "generated_at": timestamp,
    }

    if query:
        data["query"] = query
    if response_text:
        data["response"] = response_text

    if analysis_result:
        data["analysis"] = analysis_result.model_dump()

    if evidence:
        data["evidence"] = evidence.model_dump()
        data["search_metrics"] = {
            "total_results": evidence.hit_count,
            "collections_searched": evidence.total_collections_searched,
            "search_time_ms": round(evidence.search_time_ms, 1),
        }

    return json.dumps(data, indent=2, default=str)


def _create_bio_age_gauge(chronological_age: float, biological_age: float) -> Optional[bytes]:
    """Create a biological age gauge chart as PNG bytes."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import io

        import matplotlib.patches as patches
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 2.5))

        # Draw gauge background
        colors_zones = [('#22c55e', 'Younger'), ('#76B900', 'Healthy'), ('#f59e0b', 'Accelerated'), ('#ef4444', 'High Risk')]
        for i, (color, label) in enumerate(colors_zones):
            start = i * 45
            wedge = patches.Wedge((0.5, 0.0), 0.45, 180 - start - 45, 180 - start,
                                   facecolor=color, alpha=0.3, transform=ax.transAxes)
            ax.add_patch(wedge)

        # Needle position (0-180 degrees based on age acceleration)
        accel = biological_age - chronological_age
        # Map acceleration to angle: -10 -> 180 deg, 0 -> 90 deg, +10 -> 0 deg
        angle = max(0, min(180, 90 - accel * 9))
        import math
        needle_x = 0.5 + 0.35 * math.cos(math.radians(angle))
        needle_y = 0.0 + 0.35 * math.sin(math.radians(angle))
        ax.annotate('', xy=(needle_x, needle_y), xytext=(0.5, 0.0),
                     xycoords='axes fraction', textcoords='axes fraction',
                     arrowprops=dict(arrowstyle='->', color='#1f2937', lw=2))

        # Labels
        ax.text(0.5, -0.15, f'Bio Age: {biological_age:.1f}y  |  Chrono: {chronological_age:.0f}y',
                ha='center', va='center', transform=ax.transAxes, fontsize=10, fontweight='bold')
        ax.text(0.5, -0.3, f'Acceleration: {"+" if accel > 0 else ""}{accel:.1f} years',
                ha='center', va='center', transform=ax.transAxes, fontsize=9,
                color='#ef4444' if accel > 0 else '#22c55e')

        ax.set_xlim(0, 1)
        ax.set_ylim(-0.4, 0.6)
        ax.set_aspect('equal')
        ax.axis('off')

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight', transparent=True)
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()
    except Exception:
        return None


def _create_disease_risk_radar(trajectories: list) -> Optional[bytes]:
    """Create a disease risk radar chart as PNG bytes."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import io

        import matplotlib.pyplot as plt
        import numpy as np

        # Extract disease names and risk scores
        diseases = []
        scores = []
        colors = []
        risk_score_map = {"NORMAL": 10, "LOW": 25, "MODERATE": 50, "HIGH": 75, "CRITICAL": 95}
        risk_color_map = {"NORMAL": '#22c55e', "LOW": '#22c55e', "MODERATE": '#f59e0b', "HIGH": '#ef4444', "CRITICAL": '#dc2626'}

        for t in trajectories:
            if isinstance(t, dict):
                name = t.get("disease", t.get("name", "Unknown"))
                risk = str(t.get("risk_level", t.get("stage", "LOW"))).upper()
            elif hasattr(t, "disease"):
                name = t.disease
                risk = str(getattr(t, "risk_level", "LOW")).upper()
            else:
                continue
            diseases.append(name[:15])  # truncate long names
            scores.append(risk_score_map.get(risk, 25))
            colors.append(risk_color_map.get(risk, '#6b7280'))

        if not diseases:
            return None

        # Radar chart
        N = len(diseases)
        angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
        scores_plot = scores + [scores[0]]
        angles += [angles[0]]

        fig, ax = plt.subplots(figsize=(4, 4), subplot_kw=dict(polar=True))
        ax.fill(angles, scores_plot, alpha=0.2, color='#ef4444')
        ax.plot(angles, scores_plot, 'o-', color='#ef4444', linewidth=2, markersize=6)

        # Population average
        pop_avg = [20] * N + [20]
        ax.fill(angles, pop_avg, alpha=0.1, color='#76B900')
        ax.plot(angles, pop_avg, '--', color='#76B900', linewidth=1, label='Pop. Avg')

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(diseases, size=8)
        ax.set_ylim(0, 100)
        ax.set_yticks([25, 50, 75])
        ax.set_yticklabels(['Low', 'Moderate', 'High'], size=7, color='gray')
        ax.set_title('Disease Risk Profile', size=11, fontweight='bold', pad=15)
        ax.legend(loc='lower right', fontsize=7)

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()
    except Exception:
        return None


def _create_pgx_bar_chart(pgx_results: list) -> Optional[bytes]:
    """Create a PGx metabolizer status bar chart as PNG bytes."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import io

        import matplotlib.pyplot as plt

        genes = []
        phenotype_scores = []
        bar_colors = []

        phenotype_score_map = {
            "poor metabolizer": 1, "poor": 1, "poor expresser": 1,
            "intermediate metabolizer": 2, "intermediate": 2,
            "normal metabolizer": 3, "normal": 3, "extensive": 3,
            "rapid metabolizer": 4, "rapid": 4,
            "ultra-rapid metabolizer": 5, "ultra-rapid": 5, "ultrarapid metabolizer": 5,
        }
        phenotype_color_map = {
            1: '#ef4444',  # poor - red
            2: '#f59e0b',  # intermediate - amber
            3: '#22c55e',  # normal - green
            4: '#3b82f6',  # rapid - blue
            5: '#dc2626',  # ultra-rapid - dark red
        }

        for r in pgx_results:
            if isinstance(r, dict):
                gene = r.get("gene", "?")
                phenotype = str(r.get("phenotype", "unknown")).lower().strip()
            elif hasattr(r, "gene"):
                gene = r.gene
                phenotype = str(getattr(r, "phenotype", "unknown")).lower().strip()
            else:
                continue

            score = phenotype_score_map.get(phenotype, 3)
            genes.append(gene)
            phenotype_scores.append(score)
            bar_colors.append(phenotype_color_map.get(score, '#6b7280'))

        if not genes:
            return None

        fig, ax = plt.subplots(figsize=(5, max(2, len(genes) * 0.5)))
        y_pos = range(len(genes))
        ax.barh(y_pos, phenotype_scores, color=bar_colors, height=0.6, edgecolor='white')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(genes, fontsize=9)
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.set_xticklabels(['Poor', 'Intermediate', 'Normal', 'Rapid', 'Ultra-Rapid'], fontsize=7)
        ax.set_xlim(0, 5.5)
        ax.set_title('Pharmacogenomic Metabolizer Profile', fontsize=11, fontweight='bold')
        ax.axvline(x=3, color='#76B900', linestyle='--', alpha=0.5, label='Normal')

        plt.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()
    except Exception:
        return None


def export_pdf(report_markdown: str, analysis: Optional[AnalysisResult] = None) -> bytes:
    """Export a markdown report as a styled PDF with optional embedded charts.

    Uses reportlab Platypus for professional formatting with HCLS AI Factory
    branding (NVIDIA green accent, dark header bars, zebra-striped tables).
    When *analysis* is provided and matplotlib is available, embeds a
    biological-age gauge, disease-risk radar, and PGx metabolizer bar chart.
    Falls back to plain-text PDF if reportlab is not installed.

    Args:
        report_markdown: Complete markdown report string.
        analysis: Optional AnalysisResult with patient analysis data for charts.

    Returns:
        PDF content as bytes.
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            BaseDocTemplate,
            Frame,
            KeepTogether,
            PageTemplate,
            Paragraph,
            Spacer,
            Table,
            TableStyle,
        )

        # Color palette (matches drug-discovery-pipeline branding)
        NVIDIA_GREEN = colors.HexColor("#76B900")
        DARK_BG = colors.HexColor("#1B1B2F")
        HEADER_ROW_BG = colors.HexColor("#2d2d44")
        ZEBRA_LIGHT = colors.HexColor("#f5f5f5")
        BORDER_COLOR = colors.HexColor("#dddddd")
        GRID_COLOR = colors.HexColor("#eeeeee")
        WHITE = colors.white
        TEXT_PRIMARY = colors.HexColor("#1E293B")
        TEXT_MUTED = colors.HexColor("#94A3B8")

        PAGE_W, PAGE_H = letter
        MARGIN = 0.65 * inch
        CONTENT_W = PAGE_W - 2 * MARGIN

        buf = io.BytesIO()

        def _first_page(canvas, doc):
            canvas.saveState()
            canvas.setFillColor(DARK_BG)
            canvas.rect(0, PAGE_H - 62, PAGE_W, 62, fill=1, stroke=0)
            canvas.setFillColor(NVIDIA_GREEN)
            canvas.rect(0, PAGE_H - 65, PAGE_W, 3, fill=1, stroke=0)
            canvas.setFillColor(WHITE)
            canvas.setFont("Helvetica-Bold", 18)
            canvas.drawString(MARGIN, PAGE_H - 42,
                              "Biomarker Intelligence Report")
            canvas.setFont("Helvetica", 8)
            canvas.setFillColor(NVIDIA_GREEN)
            canvas.drawRightString(PAGE_W - MARGIN, PAGE_H - 30, f"v{VERSION}")
            canvas.setFillColor(NVIDIA_GREEN)
            canvas.rect(0, 30, PAGE_W, 2, fill=1, stroke=0)
            canvas.setFont("Helvetica", 7)
            canvas.setFillColor(TEXT_MUTED)
            canvas.drawString(MARGIN, 18,
                              "HCLS AI Factory -- Biomarker Intelligence Agent")
            canvas.drawRightString(PAGE_W - MARGIN, 18, f"Page {doc.page}")
            canvas.restoreState()

        def _later_pages(canvas, doc):
            canvas.saveState()
            canvas.setFillColor(DARK_BG)
            canvas.rect(0, PAGE_H - 28, PAGE_W, 28, fill=1, stroke=0)
            canvas.setFillColor(NVIDIA_GREEN)
            canvas.rect(0, PAGE_H - 30, PAGE_W, 2, fill=1, stroke=0)
            canvas.setFillColor(WHITE)
            canvas.setFont("Helvetica-Bold", 9)
            canvas.drawString(MARGIN, PAGE_H - 19, "Biomarker Intelligence Report")
            canvas.setFillColor(NVIDIA_GREEN)
            canvas.rect(0, 30, PAGE_W, 2, fill=1, stroke=0)
            canvas.setFont("Helvetica", 7)
            canvas.setFillColor(TEXT_MUTED)
            canvas.drawRightString(PAGE_W - MARGIN, 18, f"Page {doc.page}")
            canvas.restoreState()

        frame_first = Frame(MARGIN, MARGIN + 20, CONTENT_W,
                            PAGE_H - 2 * MARGIN - 65, id="first")
        frame_later = Frame(MARGIN, MARGIN + 20, CONTENT_W,
                            PAGE_H - 2 * MARGIN - 30, id="later")

        doc = BaseDocTemplate(buf, pagesize=letter)
        doc.addPageTemplates([
            PageTemplate(id="first", frames=[frame_first], onPage=_first_page),
            PageTemplate(id="later", frames=[frame_later], onPage=_later_pages),
        ])

        # -- Custom styles ------------------------------------------------
        styles = getSampleStyleSheet()

        styles.add(ParagraphStyle(
            name="SectionHeader",
            parent=styles["Heading2"],
            fontSize=14,
            spaceBefore=16,
            spaceAfter=8,
            textColor=NVIDIA_GREEN,
            fontName="Helvetica-Bold",
        ))

        styles.add(ParagraphStyle(
            name="SubSection",
            parent=styles["Heading3"],
            fontSize=11,
            spaceBefore=12,
            spaceAfter=6,
            textColor=TEXT_PRIMARY,
            fontName="Helvetica-Bold",
        ))

        styles.add(ParagraphStyle(
            name="CellText",
            parent=styles["BodyText"],
            fontSize=8,
            leading=10,
            fontName="Helvetica",
        ))

        styles.add(ParagraphStyle(
            name="CellBold",
            parent=styles["BodyText"],
            fontSize=8,
            leading=10,
            fontName="Helvetica-Bold",
        ))

        def _md_to_rl(text: str) -> str:
            """Convert markdown inline formatting to reportlab XML."""
            text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
            text = re.sub(r'\*(.+?)\*', r'<i>\1</i>', text)
            return text

        def _build_table(rows: List[List[str]]) -> Table:
            """Build a styled reportlab Table from parsed markdown rows."""
            cell_style = styles["CellText"]
            header_style = styles["CellBold"]

            # Convert cells to Paragraphs for text wrapping
            table_data = []
            for row_idx, row in enumerate(rows):
                style = header_style if row_idx == 0 else cell_style
                table_data.append([
                    Paragraph(_md_to_rl(cell.strip()), style) for cell in row
                ])

            ncols = len(rows[0]) if rows else 1
            col_width = CONTENT_W / ncols
            tbl = Table(table_data, colWidths=[col_width] * ncols,
                        repeatRows=1)

            # Styling
            cmds = [
                # Header row
                ("BACKGROUND", (0, 0), (-1, 0), HEADER_ROW_BG),
                ("TEXTCOLOR", (0, 0), (-1, 0), WHITE),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, 0), 8),
                # Data rows
                ("FONTSIZE", (0, 1), (-1, -1), 8),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                # Borders
                ("BOX", (0, 0), (-1, -1), 0.75, BORDER_COLOR),
                ("INNERGRID", (0, 0), (-1, -1), 0.5, GRID_COLOR),
                ("LINEBELOW", (0, 0), (-1, 0), 1, NVIDIA_GREEN),
            ]
            # Zebra stripe data rows
            for i in range(1, len(rows)):
                if i % 2 == 0:
                    cmds.append(("BACKGROUND", (0, i), (-1, i), ZEBRA_LIGHT))
            tbl.setStyle(TableStyle(cmds))
            return tbl

        # -- Parse markdown into PDF story ---------------------------------
        story: List = []
        lines = report_markdown.split("\n")
        i = 0
        while i < len(lines):
            stripped = lines[i].strip()

            # Blank line
            if not stripped:
                story.append(Spacer(1, 6))
                i += 1
                continue

            # H1 title
            if stripped.startswith("# ") and not stripped.startswith("## "):
                story.append(Paragraph(_md_to_rl(stripped[2:]), styles["Title"]))
                story.append(Spacer(1, 12))
                i += 1
                continue

            # H2 section header
            if stripped.startswith("## "):
                story.append(Spacer(1, 4))
                story.append(Paragraph(_md_to_rl(stripped[3:]),
                                       styles["SectionHeader"]))
                story.append(Spacer(1, 4))
                i += 1
                continue

            # H3 subsection
            if stripped.startswith("### "):
                story.append(Paragraph(_md_to_rl(stripped[4:]),
                                       styles["SubSection"]))
                story.append(Spacer(1, 3))
                i += 1
                continue

            # Markdown table: collect consecutive | rows
            if stripped.startswith("|"):
                table_rows: List[List[str]] = []
                while i < len(lines) and lines[i].strip().startswith("|"):
                    row_text = lines[i].strip()
                    # Skip separator rows (|---|---|)
                    if re.match(r'^\|[\s\-:|]+\|$', row_text):
                        i += 1
                        continue
                    cells = [c.strip() for c in re.split(r'(?<!\\)\|', row_text)]
                    # Remove empty first/last from leading/trailing |
                    if cells and cells[0] == "":
                        cells = cells[1:]
                    if cells and cells[-1] == "":
                        cells = cells[:-1]
                    if cells:
                        table_rows.append(cells)
                    i += 1
                if table_rows:
                    story.append(KeepTogether([
                        _build_table(table_rows),
                        Spacer(1, 8),
                    ]))
                continue

            # Bullet list items
            if stripped.startswith("- ") or stripped.startswith("* "):
                bullet_text = _md_to_rl(stripped[2:])
                story.append(Paragraph(f"\u2022  {bullet_text}",
                                       styles["BodyText"]))
                i += 1
                continue

            # Numbered list items
            if re.match(r'^\d+\.\s', stripped):
                text = _md_to_rl(re.sub(r'^\d+\.\s', '', stripped))
                num = stripped.split(".")[0]
                story.append(Paragraph(f"{num}.  {text}",
                                       styles["BodyText"]))
                i += 1
                continue

            # Horizontal rule
            if stripped.startswith("---"):
                story.append(Spacer(1, 10))
                i += 1
                continue

            # Regular paragraph
            text = _md_to_rl(stripped)
            story.append(Paragraph(text, styles["BodyText"]))
            i += 1

        # -- Embed charts if analysis data is available --------------------
        if analysis is not None:
            import io as _io

            from reportlab.platypus import Image as RLImage

            # Biological age gauge
            if hasattr(analysis, 'biological_age'):
                bio_age = analysis.biological_age
                if hasattr(bio_age, 'biological_age') and hasattr(bio_age, 'chronological_age'):
                    gauge_bytes = _create_bio_age_gauge(bio_age.chronological_age, bio_age.biological_age)
                elif isinstance(bio_age, dict):
                    gauge_bytes = _create_bio_age_gauge(bio_age.get('chronological_age', 0), bio_age.get('biological_age', 0))
                else:
                    gauge_bytes = None
                if gauge_bytes:
                    story.append(Spacer(1, 12))
                    story.append(RLImage(_io.BytesIO(gauge_bytes), width=350, height=175))
                    story.append(Spacer(1, 8))

            # Disease risk radar
            if hasattr(analysis, 'disease_trajectories'):
                radar_bytes = _create_disease_risk_radar(analysis.disease_trajectories)
                if radar_bytes:
                    story.append(Spacer(1, 12))
                    story.append(RLImage(_io.BytesIO(radar_bytes), width=280, height=280))
                    story.append(Spacer(1, 8))

            # PGx bar chart
            if hasattr(analysis, 'pgx_results'):
                pgx_bytes = _create_pgx_bar_chart(analysis.pgx_results)
                if pgx_bytes:
                    story.append(Spacer(1, 12))
                    story.append(RLImage(_io.BytesIO(pgx_bytes), width=350, height=max(140, len(analysis.pgx_results) * 35)))
                    story.append(Spacer(1, 8))

        doc.build(story)
        return buf.getvalue()

    except ImportError:
        # Fallback: plain text as bytes
        return report_markdown.encode("utf-8")


def export_csv(analysis: AnalysisResult) -> bytes:
    """Export analysis result as CSV for spreadsheet analysis.

    Produces labeled sections: Patient Info, Biological Age, Disease
    Trajectories, Pharmacogenomic Profile, Genotype Adjustments, and
    Critical Alerts.

    Args:
        analysis: AnalysisResult from patient analysis.

    Returns:
        CSV content as UTF-8 bytes.
    """
    buf = io.StringIO()
    writer = csv.writer(buf)
    profile = analysis.patient_profile
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # -- Patient Info --
    writer.writerow(["Section", "Field", "Value"])
    writer.writerow(["Patient Info", "Patient ID", profile.patient_id])
    writer.writerow(["Patient Info", "Age", profile.age])
    writer.writerow(["Patient Info", "Sex", profile.sex])
    writer.writerow(["Patient Info", "Generated", timestamp])
    writer.writerow([])

    # -- Biological Age --
    ba = analysis.biological_age
    writer.writerow(["Biological Age", "Chronological Age", ba.chronological_age])
    writer.writerow(["Biological Age", "Biological Age (PhenoAge)",
                     f"{ba.biological_age:.1f}"])
    writer.writerow(["Biological Age", "Age Acceleration",
                     f"{ba.age_acceleration:+.1f}"])
    writer.writerow(["Biological Age", "Mortality Risk",
                     f"{ba.mortality_risk:.4f}"])
    if ba.grimage_score is not None:
        writer.writerow(["Biological Age", "GrimAge Surrogate",
                         f"{ba.grimage_score:.1f}"])
    writer.writerow([])

    # -- Aging Drivers --
    if ba.aging_drivers:
        writer.writerow(["Aging Drivers", "Biomarker", "Value", "Direction",
                         "Contribution"])
        for d in ba.aging_drivers:
            writer.writerow([
                "",
                d.get("biomarker", d.get("marker", "")),
                d.get("value", ""),
                d.get("direction", ""),
                d.get("contribution", ""),
            ])
        writer.writerow([])

    # -- Disease Trajectories --
    writer.writerow(["Disease Trajectories", "Disease", "Risk Level",
                     "Est. Years to Onset", "Genetic Risk Factors",
                     "Interventions"])
    for traj in analysis.disease_trajectories:
        years = (f"{traj.years_to_onset_estimate:.0f}"
                 if traj.years_to_onset_estimate else "N/A")
        factors = "; ".join(traj.genetic_risk_factors) or "None"
        interventions = "; ".join(traj.intervention_recommendations) or "None"
        writer.writerow(["", traj.disease.value, traj.risk_level.value,
                         years, factors, interventions])
    writer.writerow([])

    # -- PGx Profile --
    if analysis.pgx_results:
        writer.writerow(["PGx Profile", "Gene", "Star Alleles", "Phenotype",
                         "Drugs Affected"])
        for pgx in analysis.pgx_results:
            drugs = "; ".join(
                d.get("drug", "") for d in pgx.drugs_affected
            ) or "None"
            writer.writerow(["", pgx.gene, pgx.star_alleles,
                             pgx.phenotype.value, drugs])
        writer.writerow([])

    # -- Genotype Adjustments --
    if analysis.genotype_adjustments:
        writer.writerow(["Genotype Adjustments", "Biomarker", "Gene",
                         "Genotype", "Standard Range", "Adjusted Range",
                         "Rationale"])
        for adj in analysis.genotype_adjustments:
            writer.writerow(["", adj.biomarker, adj.gene, adj.genotype,
                             adj.standard_range, adj.adjusted_range,
                             adj.rationale])
        writer.writerow([])

    # -- Critical Alerts --
    if analysis.critical_alerts:
        writer.writerow(["Critical Alerts", "Alert"])
        for alert in analysis.critical_alerts:
            writer.writerow(["", alert])
        writer.writerow([])

    # -- Biomarker Values --
    writer.writerow(["Biomarker Values", "Biomarker", "Value"])
    for name, val in sorted(profile.biomarkers.items()):
        writer.writerow(["", name, val])

    return buf.getvalue().encode("utf-8")


def validate_fhir_bundle(bundle: dict) -> List[str]:
    """Validate that a FHIR R4 Bundle conforms to DiagnosticReport structure.

    Performs lightweight structural validation without requiring external FHIR
    validator libraries.  Checks resource types, required fields, and internal
    reference integrity.

    Args:
        bundle: Parsed FHIR Bundle dictionary.

    Returns:
        List of validation error strings.  Empty list means valid.
    """
    errors: List[str] = []

    # --- Top-level Bundle checks ---
    if not isinstance(bundle, dict):
        errors.append("Bundle must be a JSON object")
        return errors

    if bundle.get("resourceType") != "Bundle":
        errors.append(
            f"resourceType must be 'Bundle', got '{bundle.get('resourceType')}'"
        )

    entries = bundle.get("entry")
    if not isinstance(entries, list):
        errors.append("Bundle must contain an 'entry' list")
        return errors  # nothing more to validate

    # Collect all fullUrls and resource IDs so we can verify references
    known_refs: set = set()
    for idx, entry in enumerate(entries):
        full_url = entry.get("fullUrl")
        if full_url:
            known_refs.add(full_url)
        resource = entry.get("resource") if isinstance(entry, dict) else None
        if resource and isinstance(resource, dict):
            res_type = resource.get("resourceType", "")
            res_id = resource.get("id", "")
            if res_type and res_id:
                known_refs.add(f"{res_type}/{res_id}")
                known_refs.add(f"urn:uuid:{res_id}")

    VALID_RESOURCE_TYPES = {
        "DiagnosticReport", "Observation", "Patient", "Practitioner",
        "Organization", "Specimen", "Condition", "MedicationStatement",
        "ServiceRequest", "Bundle", "Composition",
    }

    # Required fields per resource type
    DIAGNOSTIC_REPORT_REQUIRED = {"status", "code", "subject", "effectiveDateTime", "result"}
    OBSERVATION_REQUIRED = {"status", "code"}
    PATIENT_REQUIRED_CHECK = "identifier"

    def _check_reference(ref_value: str, location: str) -> None:
        """Verify that a reference string points to a known resource."""
        if ref_value and ref_value not in known_refs:
            errors.append(
                f"{location}: reference '{ref_value}' does not match any "
                f"fullUrl or resource id in the bundle"
            )

    def _collect_references(obj: Any, path: str) -> None:
        """Recursively find all 'reference' fields and validate them."""
        if isinstance(obj, dict):
            if "reference" in obj and isinstance(obj["reference"], str):
                _check_reference(obj["reference"], path)
            for key, value in obj.items():
                _collect_references(value, f"{path}.{key}")
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                _collect_references(item, f"{path}[{i}]")

    for idx, entry in enumerate(entries):
        prefix = f"entry[{idx}]"

        if not isinstance(entry, dict):
            errors.append(f"{prefix}: entry must be a JSON object")
            continue

        resource = entry.get("resource")
        if not isinstance(resource, dict):
            errors.append(f"{prefix}: entry must contain a 'resource' object")
            continue

        res_type = resource.get("resourceType")
        if not res_type:
            errors.append(f"{prefix}: resource is missing 'resourceType'")
            continue

        if res_type not in VALID_RESOURCE_TYPES:
            errors.append(
                f"{prefix}: unknown resourceType '{res_type}'"
            )

        # --- DiagnosticReport checks ---
        if res_type == "DiagnosticReport":
            for field in DIAGNOSTIC_REPORT_REQUIRED:
                if field not in resource:
                    errors.append(
                        f"{prefix} (DiagnosticReport): missing required field '{field}'"
                    )
            # subject must have reference
            subject = resource.get("subject")
            if isinstance(subject, dict) and "reference" not in subject:
                errors.append(
                    f"{prefix} (DiagnosticReport): 'subject' must contain a 'reference'"
                )
            # result must be a list
            result = resource.get("result")
            if result is not None and not isinstance(result, list):
                errors.append(
                    f"{prefix} (DiagnosticReport): 'result' must be a list"
                )

        # --- Observation checks ---
        elif res_type == "Observation":
            for field in OBSERVATION_REQUIRED:
                if field not in resource:
                    errors.append(
                        f"{prefix} (Observation): missing required field '{field}'"
                    )
            # Must have valueQuantity or valueString
            has_value = (
                "valueQuantity" in resource
                or "valueString" in resource
            )
            if not has_value:
                errors.append(
                    f"{prefix} (Observation): must have 'valueQuantity' or 'valueString'"
                )

        # --- Patient checks ---
        elif res_type == "Patient":
            if PATIENT_REQUIRED_CHECK not in resource:
                errors.append(
                    f"{prefix} (Patient): missing required field 'identifier'"
                )

        # --- Reference integrity for this entry ---
        _collect_references(resource, prefix)

    return errors


def export_fhir_diagnostic_report(
    analysis: AnalysisResult,
    patient_profile: PatientProfile,
    validate: bool = True,
) -> str:
    """Export analysis result as a FHIR R4 DiagnosticReport JSON bundle.

    Creates a FHIR Bundle containing:
    - DiagnosticReport resource (main report)
    - Observation resources (biological age, disease trajectories, PGx)
    - Patient resource reference

    Args:
        analysis: AnalysisResult from patient analysis.
        patient_profile: PatientProfile with patient demographics.
        validate: If True (default), run lightweight FHIR R4 structural
            validation on the generated bundle and log any warnings.
            Validation issues are logged but never cause the export to fail.

    Returns:
        FHIR R4 JSON string (Bundle).
    """
    timestamp = datetime.now(timezone.utc).isoformat()
    report_id = f"biomarker-report-{patient_profile.patient_id}"

    entries = []

    # Patient resource reference
    patient_ref = f"Patient/{patient_profile.patient_id}"

    # DiagnosticReport resource
    diagnostic_report = {
        "resource": {
            "resourceType": "DiagnosticReport",
            "id": report_id,
            "status": "final",
            "category": [
                {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                            "code": "GE",
                            "display": "Genetics",
                        }
                    ]
                }
            ],
            "code": {
                "coding": [
                    {
                        "system": "http://loinc.org",
                        "code": "51969-4",
                        "display": "Genetic analysis report",
                    }
                ],
                "text": "Biomarker Intelligence Report",
            },
            "subject": {"reference": patient_ref},
            "effectiveDateTime": timestamp,
            "issued": timestamp,
            "performer": [
                {
                    "display": "HCLS AI Factory - Biomarker Intelligence Agent",
                }
            ],
            "result": [],
            "conclusion": "",
        },
        "fullUrl": f"urn:uuid:{report_id}",
    }

    observation_refs = []

    # Observation: Biological Age
    bio_age_id = f"observation-bio-age-{patient_profile.patient_id}"
    bio_age_obs = {
        "resource": {
            "resourceType": "Observation",
            "id": bio_age_id,
            "status": "final",
            "category": [
                {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Laboratory",
                        }
                    ]
                }
            ],
            "code": {
                "coding": [
                    {
                        "system": "http://loinc.org",
                        "code": "88331-4",
                        "display": "Biological age",
                    }
                ],
                "text": "Biological Age (PhenoAge)",
            },
            "subject": {"reference": patient_ref},
            "effectiveDateTime": timestamp,
            "valueQuantity": {
                "value": analysis.biological_age.biological_age,
                "unit": "years",
                "system": "http://unitsofmeasure.org",
                "code": "a",
            },
            "component": [
                {
                    "code": {"text": "Chronological Age"},
                    "valueQuantity": {
                        "value": analysis.biological_age.chronological_age,
                        "unit": "years",
                        "system": "http://unitsofmeasure.org",
                        "code": "a",
                    },
                },
                {
                    "code": {"text": "Age Acceleration"},
                    "valueQuantity": {
                        "value": analysis.biological_age.age_acceleration,
                        "unit": "years",
                        "system": "http://unitsofmeasure.org",
                        "code": "a",
                    },
                },
                {
                    "code": {"text": "Mortality Risk Score"},
                    "valueQuantity": {
                        "value": analysis.biological_age.mortality_risk,
                        "unit": "{score}",
                        "system": "http://unitsofmeasure.org",
                        "code": "{score}",
                    },
                },
            ],
        },
        "fullUrl": f"urn:uuid:{bio_age_id}",
    }
    entries.append(bio_age_obs)
    observation_refs.append({"reference": f"urn:uuid:{bio_age_id}"})

    # Observations: Disease Trajectories
    for i, traj in enumerate(analysis.disease_trajectories):
        traj_id = f"observation-trajectory-{traj.disease.value}-{patient_profile.patient_id}"
        traj_obs = {
            "resource": {
                "resourceType": "Observation",
                "id": traj_id,
                "status": "final",
                "category": [
                    {
                        "coding": [
                            {
                                "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                "code": "exam",
                                "display": "Exam",
                            }
                        ]
                    }
                ],
                "code": {
                    "coding": [
                        {
                            "system": "http://snomed.info/sct",
                            "code": "225338004",
                            "display": "Risk assessment",
                        }
                    ],
                    "text": f"Disease Risk Assessment - {traj.disease.value}",
                },
                "subject": {"reference": patient_ref},
                "effectiveDateTime": timestamp,
                "valueString": traj.risk_level.value,
                "interpretation": [
                    {
                        "coding": [
                            {
                                "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                                "code": "H" if traj.risk_level.value in ("high", "critical") else "N",
                                "display": "High" if traj.risk_level.value in ("high", "critical") else "Normal",
                            }
                        ]
                    }
                ],
            },
            "fullUrl": f"urn:uuid:{traj_id}",
        }
        if traj.years_to_onset_estimate:
            traj_obs["resource"]["component"] = [
                {
                    "code": {"text": "Estimated Years to Clinical Onset"},
                    "valueQuantity": {
                        "value": traj.years_to_onset_estimate,
                        "unit": "years",
                    },
                }
            ]
        entries.append(traj_obs)
        observation_refs.append({"reference": f"urn:uuid:{traj_id}"})

    # Observations: PGx Results
    for pgx in analysis.pgx_results:
        pgx_id = f"observation-pgx-{pgx.gene}-{patient_profile.patient_id}"
        pgx_obs = {
            "resource": {
                "resourceType": "Observation",
                "id": pgx_id,
                "status": "final",
                "category": [
                    {
                        "coding": [
                            {
                                "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                "code": "laboratory",
                                "display": "Laboratory",
                            }
                        ]
                    }
                ],
                "code": {
                    "coding": [
                        {
                            "system": "http://loinc.org",
                            "code": "51963-7",
                            "display": "Medication assessment",
                        }
                    ],
                    "text": f"Pharmacogenomic Result - {pgx.gene}",
                },
                "subject": {"reference": patient_ref},
                "effectiveDateTime": timestamp,
                "valueString": f"{pgx.star_alleles} ({pgx.phenotype.value})",
                "component": [
                    {
                        "code": {"text": "Gene"},
                        "valueString": pgx.gene,
                    },
                    {
                        "code": {"text": "Star Alleles"},
                        "valueString": pgx.star_alleles,
                    },
                    {
                        "code": {"text": "Metabolizer Phenotype"},
                        "valueString": pgx.phenotype.value,
                    },
                ],
            },
            "fullUrl": f"urn:uuid:{pgx_id}",
        }

        # Add drug-specific recommendations as components
        for drug_info in pgx.drugs_affected[:5]:
            pgx_obs["resource"]["component"].append({
                "code": {"text": f"Drug Recommendation - {drug_info.get('drug', '')}"},
                "valueString": drug_info.get("recommendation", ""),
            })

        entries.append(pgx_obs)
        observation_refs.append({"reference": f"urn:uuid:{pgx_id}"})

    # Set result references and conclusion on DiagnosticReport
    diagnostic_report["resource"]["result"] = observation_refs

    # Build conclusion from critical alerts
    conclusion_parts = []
    if analysis.critical_alerts:
        conclusion_parts.extend(analysis.critical_alerts)
    conclusion_parts.append(
        f"Biological age: {analysis.biological_age.biological_age:.1f} years "
        f"(acceleration: {analysis.biological_age.age_acceleration:+.1f} years)."
    )
    diagnostic_report["resource"]["conclusion"] = " | ".join(conclusion_parts)

    # Patient resource (required for reference integrity)
    patient_resource = {
        "resource": {
            "resourceType": "Patient",
            "id": patient_profile.patient_id,
            "identifier": [
                {
                    "system": "urn:hcls-ai-factory:patient-id",
                    "value": patient_profile.patient_id,
                }
            ],
            "gender": "male" if patient_profile.sex == "M" else "female",
            "birthDate": str(
                datetime.now().year - patient_profile.age
            ) if patient_profile.age else None,
        },
        "fullUrl": f"urn:uuid:{patient_profile.patient_id}",
    }
    entries.append(patient_resource)

    # Insert DiagnosticReport as first entry
    entries.insert(0, diagnostic_report)

    # Build FHIR Bundle
    bundle = {
        "resourceType": "Bundle",
        "id": f"bundle-{report_id}",
        "meta": {
            "lastUpdated": timestamp,
        },
        "type": "collection",
        "timestamp": timestamp,
        "entry": entries,
    }

    # Optional FHIR R4 structural validation
    if validate:
        validation_errors = validate_fhir_bundle(bundle)
        if validation_errors:
            logger.warning(
                "FHIR R4 validation found {} issue(s) in DiagnosticReport bundle:",
                len(validation_errors),
            )
            for err in validation_errors:
                logger.warning("  FHIR validation: {}", err)

    return json.dumps(bundle, indent=2, default=str)


# =====================================================================
# PRIVATE HELPERS
# =====================================================================


def _format_filters(filters_applied: Optional[dict]) -> str:
    """Format sidebar filters for display."""
    if not filters_applied:
        return "None"
    parts = []
    for key, value in filters_applied.items():
        if value:
            parts.append(f"{key}: {value}")
    return ", ".join(parts) if parts else "None"


def _format_citation_link(collection: str, record_id: str) -> str:
    """Format a clickable citation link.

    Note: citation formatting logic shared with rag_engine._format_citation()
    """
    if collection == "ClinicalEvidence" and record_id.isdigit():
        return f"[PMID {record_id}](https://pubmed.ncbi.nlm.nih.gov/{record_id}/)"
    return record_id


def _format_evidence_section(evidence: CrossCollectionResult) -> List[str]:
    """Format all evidence grouped by collection."""
    lines = []
    by_coll = evidence.hits_by_collection()
    for coll_name, hits in by_coll.items():
        lines.extend(_format_evidence_table(hits, coll_name))
        lines.append("")
    return lines


def _format_evidence_table(hits: List[SearchHit], collection_name: str) -> List[str]:
    """Format a Markdown table for hits from a single collection."""
    lines = [f"### {collection_name} ({len(hits)} results)", ""]

    if collection_name in ("BiomarkerRef", "AgingMarker", "GenotypeAdj", "Monitoring"):
        lines.append("| # | ID | Score | Text |")
        lines.append("|---|-----|-------|------|")
        for i, hit in enumerate(hits[:10], 1):
            text = hit.text[:100].replace("|", "\\|")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {text} |")

    elif collection_name == "ClinicalEvidence":
        lines.append("| # | ID | Score | Source | Title | Year | Disease Area |")
        lines.append("|---|-----|-------|--------|-------|------|-------------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            link = _format_citation_link(hit.collection, hit.id)
            title = m.get("title", "")[:60]
            year = m.get("year", "")
            area = m.get("disease_area", "")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {link} | {title} | {year} | {area} |")

    elif collection_name in ("GeneticVariant", "PGxRule", "DrugInteraction"):
        lines.append("| # | ID | Score | Gene | Drug/Variant | Text |")
        lines.append("|---|-----|-------|------|-------------|------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            gene = m.get("gene", "")
            drug_var = m.get("drug", m.get("rs_id", ""))
            text = hit.text[:80].replace("|", "\\|")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {gene} | {drug_var} | {text} |")

    elif collection_name == "DiseaseTrajectory":
        lines.append("| # | ID | Score | Disease | Stage | Text |")
        lines.append("|---|-----|-------|---------|-------|------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            disease = m.get("disease", "")
            stage = m.get("stage", "")
            text = hit.text[:80].replace("|", "\\|")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {disease} | {stage} | {text} |")

    elif collection_name == "Nutrition":
        lines.append("| # | ID | Score | Nutrient | Genetic Context | Text |")
        lines.append("|---|-----|-------|----------|-----------------|------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            nutrient = m.get("nutrient", "")
            genetic = m.get("genetic_context", "")[:30]
            text = hit.text[:80].replace("|", "\\|")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {nutrient} | {genetic} | {text} |")

    elif collection_name == "Genomic":
        lines.append("| # | ID | Score | Gene | Consequence | Clinical Significance |")
        lines.append("|---|-----|-------|------|-------------|----------------------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            gene = m.get("gene", "")
            consequence = m.get("consequence", "")[:25]
            clin_sig = m.get("clinical_significance", "")[:25]
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {gene} | {consequence} | {clin_sig} |")

    else:
        # Generic fallback
        lines.append("| # | ID | Score | Text |")
        lines.append("|---|-----|-------|------|")
        for i, hit in enumerate(hits[:10], 1):
            text = hit.text[:100].replace("|", "\\|")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {text} |")

    return lines
