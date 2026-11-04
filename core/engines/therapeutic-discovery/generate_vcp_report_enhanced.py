#!/usr/bin/env python3
"""
VCP Drug Candidate Report Generator - Enhanced Edition

Renders a pipeline run as a PDF. It renders ONLY what it is given.

Until 2026-11-04 this module was a hard-coded result set wearing the names of real tools: the
"Top Ranked Drug Candidates" table held four literal dicts (VCP-AI-001..004, docking -8.62 to
-10.95 kcal/mol, QED, MW, logP, composite score) under the prose "NVIDIA BioNeMo MolMIM generated
novel molecules ... docked against VCP using DiffDock"; the variant section attributed an invented
variant plus AlphaMissense/CADD/ClinVar calls to "NVIDIA Parabricks 4.6 ... identifying"; the
structure section stated invented binding-site evidence (key residues, pocket volume ~450 A^3,
druggability 0.92); the executive summary quoted those docking scores back and printed a
"Pipeline Performance" timing table nobody had measured. The file imports no NIM client, no
docking code and no RDKit descriptor; nothing in it was ever computed. It is reachable from two
UI buttons (:8505 Drug Discovery, :8501 Chat), so the PDF reads as a pipeline result.

Every one of those literals is gone. The sections that report computed results now take their
numbers from an injected `results` dict (the shape `run_cloud_nim_report.py` already writes) and
raise `ReportDataUnavailable` when it is absent, rather than falling back to invention.
"""

import json
import os
import requests
import hashlib
from datetime import datetime
from pathlib import Path
from io import BytesIO

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, Image, HRFlowable, KeepTogether, Flowable
)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY, TA_RIGHT
from reportlab.graphics.shapes import Drawing, Rect, String, Line
from reportlab.graphics import renderPDF


# =============================================================================
# Color Scheme (matching HCLS AI Factory landing page)
# =============================================================================
class Colors:
    """Dark theme color palette."""
    # Primary
    NVIDIA_GREEN = colors.HexColor('#76B900')
    NVIDIA_GREEN_LIGHT = colors.HexColor('#8DD100')
    NVIDIA_GREEN_DARK = colors.HexColor('#5A8F00')

    # Dark backgrounds
    BG_PRIMARY = colors.HexColor('#0a0a0f')
    BG_SECONDARY = colors.HexColor('#12121a')
    BG_CARD = colors.HexColor('#1a1a24')
    BG_CARD_HOVER = colors.HexColor('#222230')

    # Text
    TEXT_PRIMARY = colors.HexColor('#ffffff')
    TEXT_SECONDARY = colors.HexColor('#a0a0b0')
    TEXT_MUTED = colors.HexColor('#6a6a7a')

    # Accents
    ACCENT_BLUE = colors.HexColor('#3b82f6')
    ACCENT_PURPLE = colors.HexColor('#8b5cf6')
    ACCENT_CYAN = colors.HexColor('#06b6d4')
    ACCENT_PINK = colors.HexColor('#ec4899')
    ACCENT_ORANGE = colors.HexColor('#f97316')

    # Status
    STATUS_SUCCESS = colors.HexColor('#22c55e')
    STATUS_WARNING = colors.HexColor('#eab308')
    STATUS_ERROR = colors.HexColor('#ef4444')

    # For PDF (light background version - PDFs read better with light bg)
    PDF_BG = colors.HexColor('#0f0f17')
    PDF_CARD_BG = colors.HexColor('#1a1a24')
    PDF_HEADER_BG = colors.HexColor('#76B900')


class ReportDataUnavailable(RuntimeError):
    """A report section was asked to report results that were never computed.

    Raised instead of emitting a plausible-looking placeholder, because this PDF is read as a
    pipeline result. Feed the generator a `results` dict (see `run_cloud_nim_report.py`).
    """


class GradientRect(Flowable):
    """Custom flowable for gradient header bars."""

    def __init__(self, width, height, color1, color2):
        Flowable.__init__(self)
        self.width = width
        self.height = height
        self.color1 = color1
        self.color2 = color2

    def draw(self):
        self.canv.setFillColor(self.color1)
        self.canv.rect(0, 0, self.width, self.height, fill=1, stroke=0)


class VCPReportGeneratorEnhanced:
    """Generate stunning VCP Drug Candidate PDF Report."""

    def __init__(self, output_path: str = "outputs/VCP_Drug_Candidate_Report.pdf", context: dict = None,
                 results: dict = None):
        """`results` carries the run that is being reported. Without it the sections that would
        have to assert computed numbers raise `ReportDataUnavailable` instead of inventing them."""
        self.output_path = Path(output_path)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache_dir = Path("data/structures/image_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.context = context  # Dynamic data from latest RAG query
        self.results = results or None
        self.styles = getSampleStyleSheet()
        self._setup_styles()
        self._load_structure_data()

    # -- data access -------------------------------------------------------------------
    def _pipeline_meta(self) -> dict:
        """The run's own metadata, or {} when no run was supplied."""
        if isinstance(self.results, dict):
            meta = self.results.get("pipeline")
            if isinstance(meta, dict):
                return meta
        return {}

    def _candidates(self) -> list:
        """Ranked candidates from the injected run. Empty when none was supplied."""
        if isinstance(self.results, dict):
            cands = self.results.get("top_candidates") or self.results.get("all_candidates")
            if isinstance(cands, list):
                return cands
        return []

    def _require_results(self, section: str) -> None:
        if not self._candidates():
            raise ReportDataUnavailable(
                f"{section} reports computed results, and no run was supplied to this report "
                "generator. Pass results=<dict with 'top_candidates'> (the shape written by "
                "run_cloud_nim_report.py / the Drug Discovery UI), or run "
                "`python run_cloud_nim_report.py` to produce a real run first. This generator "
                "will not emit candidate numbers it did not receive."
            )

    def _unavailable(self, heading: str, why: str) -> list:
        """A section that cannot be filled is stated as missing, never approximated."""
        return [
            Paragraph(heading, self.styles['SectionHeader']),
            Paragraph(
                f"<b>NOT INCLUDED IN THIS REPORT.</b> {why} No value in this section was "
                "computed, estimated or inferred.",
                self.styles['BodyJustified'],
            ),
            Spacer(1, 10),
        ]

    def _load_structure_data(self):
        """Load VCP structure data."""
        struct_file = Path("data/structures/vcp_structures.json")
        if struct_file.exists():
            with open(struct_file) as f:
                self.structures = json.load(f)
        else:
            self.structures = []

    def _setup_styles(self):
        """Set up custom paragraph styles with dark theme."""
        # Main title - deep blue
        self.styles.add(ParagraphStyle(
            name='TitleMain',
            parent=self.styles['Heading1'],
            fontSize=28,
            spaceAfter=4,
            alignment=TA_CENTER,
            textColor=colors.HexColor('#0B2375'),
            fontName='Helvetica-Bold'
        ))

        self.styles.add(ParagraphStyle(
            name='TitleSub',
            parent=self.styles['Heading1'],
            fontSize=16,
            spaceAfter=20,
            alignment=TA_CENTER,
            textColor=Colors.TEXT_SECONDARY,
            fontName='Helvetica'
        ))

        self.styles.add(ParagraphStyle(
            name='SectionHeader',
            parent=self.styles['Heading2'],
            fontSize=16,
            spaceBefore=25,
            spaceAfter=12,
            textColor=Colors.NVIDIA_GREEN,
            fontName='Helvetica-Bold',
            leftIndent=0,
        ))

        self.styles.add(ParagraphStyle(
            name='SubSection',
            parent=self.styles['Heading3'],
            fontSize=12,
            spaceBefore=15,
            spaceAfter=8,
            textColor=Colors.ACCENT_CYAN,
            fontName='Helvetica-Bold'
        ))

        self.styles.add(ParagraphStyle(
            name='BodyJustified',
            parent=self.styles['Normal'],
            fontSize=10,
            spaceAfter=8,
            alignment=TA_JUSTIFY,
            leading=14,
            textColor=colors.HexColor('#333333')
        ))

        self.styles.add(ParagraphStyle(
            name='SmallText',
            parent=self.styles['Normal'],
            fontSize=8,
            textColor=Colors.TEXT_MUTED
        ))

        self.styles.add(ParagraphStyle(
            name='CodeText',
            parent=self.styles['Normal'],
            fontSize=7,
            fontName='Courier',
            textColor=colors.HexColor('#444444'),
            backColor=colors.HexColor('#f5f5f5'),
            leftIndent=5,
            rightIndent=5,
            spaceBefore=4,
            spaceAfter=4
        ))

        self.styles.add(ParagraphStyle(
            name='Caption',
            parent=self.styles['Normal'],
            fontSize=9,
            alignment=TA_CENTER,
            textColor=Colors.TEXT_MUTED,
            spaceAfter=15
        ))

        self.styles.add(ParagraphStyle(
            name='Highlight',
            parent=self.styles['Normal'],
            fontSize=11,
            textColor=Colors.NVIDIA_GREEN,
            fontName='Helvetica-Bold'
        ))

    def _download_image(self, url: str, filename: str) -> Path:
        """Download and cache an image."""
        cache_path = self.cache_dir / filename

        if cache_path.exists():
            return cache_path

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            with open(cache_path, 'wb') as f:
                f.write(response.content)
            return cache_path
        except Exception as e:
            print(f"Warning: Could not download image {url}: {e}")
            return None

    def _get_structure_image(self, pdb_id: str) -> Path:
        """Get structure image from RCSB PDB."""
        pdb_lower = pdb_id.lower()
        url = f"https://cdn.rcsb.org/images/structures/{pdb_lower[1:3]}/{pdb_lower}/{pdb_lower}_assembly-1.jpeg"
        return self._download_image(url, f"{pdb_lower}_structure.jpeg")

    def _get_molecule_image_url(self, smiles: str) -> str:
        """Get molecule image URL from PubChem."""
        # Use PubChem's structure image service
        encoded_smiles = requests.utils.quote(smiles)
        return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/PNG?image_size=300x300"

    def _download_molecule_image(self, smiles: str, name: str) -> Path:
        """Download molecule image from PubChem."""
        url = self._get_molecule_image_url(smiles)
        safe_name = "".join(c for c in name if c.isalnum() or c in "._-")
        return self._download_image(url, f"mol_{safe_name}.png")

    def _create_header(self):
        """Create stunning header section."""
        elements = []

        # Green accent bar
        elements.append(GradientRect(7*inch, 4, Colors.NVIDIA_GREEN, Colors.NVIDIA_GREEN_LIGHT))
        elements.append(Spacer(1, 20))

        # Main title (split across two lines)
        elements.append(Paragraph(
            "PRECISION MEDICINE TO DRUG",
            self.styles['TitleMain']
        ))
        elements.append(Spacer(1, 2))
        elements.append(Paragraph(
            "DISCOVERY",
            self.styles['TitleMain']
        ))
        elements.append(Spacer(1, 12))
        elements.append(Paragraph(
            "AI Factory Pipeline Report",
            self.styles['TitleSub']
        ))

        # Info box
        elements.append(Spacer(1, 10))

        # Target and sample come from the run; "VCP" and "HG002" used to be printed here for
        # every report, including reports exported for a different gene entirely.
        meta = self._pipeline_meta()
        gene = meta.get("target_gene") or (self.context or {}).get("target_gene") or "not supplied"
        sample = meta.get("sample_id") or (self.context or {}).get("sample_id") or "not supplied"

        info_data = [[
            Paragraph(f"<b>Target:</b> {self._safe(str(gene))}", self.styles['BodyJustified']),
            Paragraph(f"<b>Sample:</b> {self._safe(str(sample))}", self.styles['BodyJustified']),
            Paragraph(f"<b>Generated:</b> {datetime.now().strftime('%B %d, %Y')}", self.styles['BodyJustified']),
        ]]

        info_table = Table(info_data, colWidths=[2*inch, 2*inch, 2.5*inch])
        info_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#f0f8e8')),
            ('BOX', (0, 0), (-1, -1), 2, Colors.NVIDIA_GREEN),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('TOPPADDING', (0, 0), (-1, -1), 12),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
        ]))
        elements.append(info_table)

        return elements

    def _create_pipeline_flow(self):
        """Create visual pipeline flow."""
        elements = []

        elements.append(Spacer(1, 20))

        # Pipeline phases
        flow_data = [
            ['PHASE 1-3', 'PHASE 4', 'PHASE 5', 'PHASE 6'],
            ['GENOMICS', 'RAG/CHAT', 'STRUCTURE', 'MOLECULES'],
            ['VCP Variant\nDetected', 'Target\nValidated', 'Cryo-EM\nEvidence', 'Drug\nCandidates'],
        ]

        table = Table(flow_data, colWidths=[1.55*inch]*4)
        table.setStyle(TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            # Phase labels
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, 0), 8),
            ('TEXTCOLOR', (0, 0), (-1, 0), Colors.NVIDIA_GREEN),
            # Main labels
            ('FONTNAME', (0, 1), (-1, 1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 1), (-1, 1), 12),
            ('BACKGROUND', (0, 1), (-1, 1), Colors.NVIDIA_GREEN),
            ('TEXTCOLOR', (0, 1), (-1, 1), colors.white),
            # Descriptions
            ('FONTSIZE', (0, 2), (-1, 2), 9),
            ('TEXTCOLOR', (0, 2), (-1, 2), colors.HexColor('#444444')),
            ('BACKGROUND', (0, 2), (-1, 2), colors.HexColor('#f8fff0')),
            # Borders
            ('BOX', (0, 0), (-1, -1), 2, Colors.NVIDIA_GREEN),
            ('INNERGRID', (0, 0), (-1, -1), 1, colors.HexColor('#c8e8a0')),
            ('TOPPADDING', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 10),
        ]))
        elements.append(table)

        return elements

    def _create_variant_section(self):
        """Create genomic variant section from the supplied variant call, or state that there is none.

        This section used to hard-code a VCP variant (chr9:35,065,263 G>A, rs188935092,
        heterozygous) and a pathogenicity panel (AlphaMissense 0.89, ClinVar PATHOGENIC, CADD
        28.5) under the claim that "NVIDIA Parabricks 4.6 on DGX Spark processed the HG002 whole
        genome sample ... identifying" it. No aligner, caller or annotator runs in this file, and
        the same literals disagree with the two other places in the repo that hard-code the same
        demo variant. Nothing is asserted here that was not handed in.
        """
        variant = None
        if isinstance(self.results, dict):
            variant = self.results.get("variant")
        if not isinstance(variant, dict):
            return self._unavailable(
                "1. GENOMIC VARIANT DETECTION",
                "No variant call set was supplied to this report generator, so no variant, "
                "pathogenicity score or ClinVar classification is shown. This report ran no "
                "aligner, variant caller or annotator.",
            )

        elements = []

        elements.append(Paragraph("1. GENOMIC VARIANT DETECTION", self.styles['SectionHeader']))

        caller = variant.get("called_by")
        elements.append(Paragraph(
            "The variant below is the one supplied to this report"
            + (f", called by {self._safe(str(caller))}." if caller
               else ". The caller that produced it was not recorded, so no tool is credited for it."),
            self.styles['BodyJustified']
        ))

        # Variant details — only fields the caller actually provided
        elements.append(Paragraph("Supplied Variant", self.styles['SubSection']))

        def _v(key):
            val = variant.get(key)
            return self._safe(str(val)) if val not in (None, "") else 'not supplied'

        variant_data = [
            ['Property', 'Value', 'Property', 'Value'],
            ['Gene', _v('gene'), 'rsID', _v('rsid')],
            ['Chromosome', _v('chrom'), 'Consequence', _v('consequence')],
            ['Position', _v('pos'), 'Impact', _v('impact')],
            ['Change', f"{_v('ref')} → {_v('alt')}", 'Zygosity', _v('zygosity')],
        ]

        table = Table(variant_data, colWidths=[1.2*inch, 1.8*inch, 1.2*inch, 1.8*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), Colors.BG_CARD),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
            ('FONTNAME', (2, 1), (2, -1), 'Helvetica-Bold'),
            ('TEXTCOLOR', (0, 1), (0, -1), Colors.ACCENT_CYAN),
            ('TEXTCOLOR', (2, 1), (2, -1), Colors.ACCENT_CYAN),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f8f8')]),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor('#dddddd')),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#eeeeee')),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('LEFTPADDING', (0, 0), (-1, -1), 10),
        ]))
        elements.append(table)
        elements.append(Spacer(1, 15))

        # Pathogenicity — rendered only from annotations that came in with the variant.
        annotations = variant.get("pathogenicity")
        if not isinstance(annotations, list) or not annotations:
            elements.append(Paragraph(
                "<b>Pathogenicity:</b> not supplied. No AlphaMissense, ClinVar or CADD "
                "annotation accompanied this variant, and none was estimated.",
                self.styles['BodyJustified']
            ))
            return elements

        path_header = Paragraph("Pathogenicity Assessment (as supplied)", self.styles['SubSection'])

        path_data = [['Source', 'Score', 'Classification']]
        for ann in annotations:
            if not isinstance(ann, dict):
                continue
            path_data.append([
                self._safe(str(ann.get('source', 'not supplied'))),
                self._safe(str(ann.get('score', '—'))),
                self._safe(str(ann.get('classification', 'not supplied'))),
            ])

        path_table = Table(path_data, colWidths=[2*inch, 1.5*inch, 2.5*inch])
        path_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), Colors.STATUS_ERROR),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (2, 1), (2, -1), 'Helvetica-Bold'),
            ('TEXTCOLOR', (2, 1), (2, -1), Colors.STATUS_ERROR),
            ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#fff8f8')),
            ('BOX', (0, 0), (-1, -1), 1, Colors.STATUS_ERROR),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#ffcccc')),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ]))
        # Keep pathogenicity section together
        elements.append(KeepTogether([path_header, path_table]))

        return elements

    def _create_target_section(self):
        """Create target hypothesis section."""
        elements = []

        elements.append(Paragraph("2. TARGET BACKGROUND (PUBLISHED LITERATURE)", self.styles['SectionHeader']))

        # This paragraph used to credit a retrieval run over "3.5 million genomic evidence
        # embeddings" that this file neither performs nor measures (the platform's Milvus holds
        # four orders of magnitude fewer vectors), and the table below carried an invented
        # five-star "Priority Score". What is left is published fact about VCP/p97, labelled as
        # background rather than as this pipeline's finding.
        elements.append(Paragraph(
            "The following is published background on VCP/p97, included for context. It is not "
            "a retrieval result and not an output of this run: no embedding search, vector "
            "database or model inference contributed to this section.",
            self.styles['BodyJustified']
        ))

        # Target profile - keep together
        target_header = Paragraph("Target Profile (literature)", self.styles['SubSection'])

        target_data = [
            ['Property', 'Value'],
            ['Target Gene', 'VCP (Valosin-Containing Protein)'],
            ['Protein', 'p97 AAA+ ATPase'],
            ['UniProt', 'P55072'],
            ['Therapeutic Area', 'Neurodegeneration'],
            ['Known inhibitor chemotype', 'ATP-competitive, D2 ATPase domain (e.g. CB-5083)'],
        ]

        target_table = Table(target_data, colWidths=[2*inch, 4*inch])
        target_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), Colors.BG_CARD),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
            ('TEXTCOLOR', (0, 1), (0, -1), Colors.ACCENT_CYAN),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f8f8')]),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor('#dddddd')),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#eeeeee')),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('LEFTPADDING', (0, 0), (-1, -1), 10),
        ]))
        elements.append(KeepTogether([target_header, target_table]))
        elements.append(Spacer(1, 15))

        # Disease associations - keep together
        disease_header = Paragraph("Disease Associations (literature)", self.styles['SubSection'])

        # The third column used to be a filled-dot "Evidence" grade (●●●●○ / ●●●○○). No evidence
        # grading is computed anywhere in this repo; the dots were decoration. Column removed.
        disease_data = [
            ['Disease', 'Mechanism'],
            ['Frontotemporal Dementia (FTD)', 'Proteostasis disruption'],
            ['Amyotrophic Lateral Sclerosis', 'Motor neuron aggregates'],
            ['Inclusion Body Myopathy', 'Muscle protein QC failure'],
        ]

        disease_table = Table(disease_data, colWidths=[2.8*inch, 3.2*inch])
        disease_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), Colors.ACCENT_PURPLE),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (2, 0), (2, -1), 'CENTER'),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f8fc')]),
            ('BOX', (0, 0), (-1, -1), 1, Colors.ACCENT_PURPLE),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e0e0f0')),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
        ]))
        elements.append(KeepTogether([disease_header, disease_table]))

        return elements

    def _create_rag_query_section(self):
        """Create dynamic section showing the latest RAG chat query and evidence."""
        if not self.context:
            return []

        elements = []

        elements.append(Spacer(1, 10))
        elements.append(GradientRect(6.5*inch, 2, Colors.ACCENT_BLUE, Colors.ACCENT_CYAN))
        elements.append(Spacer(1, 8))

        timestamp = self.context.get('timestamp', 'Unknown')
        model = self.context.get('model', 'Unknown')
        total_queries = self.context.get('total_queries', 0)

        elements.append(Paragraph(
            f"LIVE RAG QUERY — {timestamp}",
            self.styles['SectionHeader']
        ))

        elements.append(Paragraph(
            f"This section reflects the most recent query processed by the RAG/Chat pipeline. "
            f"Model: {self._safe(model)} | Session queries: {total_queries}",
            self.styles['BodyJustified']
        ))

        # User query
        query = self.context.get('query', '')
        if query:
            query_header = Paragraph("User Query", self.styles['SubSection'])
            query_text = Paragraph(
                f"<i>\"{self._safe(query)}\"</i>",
                self.styles['BodyJustified']
            )
            elements.append(KeepTogether([query_header, query_text]))
            elements.append(Spacer(1, 10))

        # Evidence table
        evidence = self.context.get('evidence', [])
        if evidence:
            ev_header = Paragraph(f"Retrieved Evidence ({len(evidence)} variants)", self.styles['SubSection'])

            ev_data = [['Gene', 'rsID', 'Consequence', 'Impact', 'Score', 'AlphaMissense']]
            for ev in evidence[:10]:  # Cap at 10 rows
                score = ev.get('score', 0)
                score_str = f"{score:.2f}" if isinstance(score, (int, float)) else str(score)
                am = ev.get('am_pathogenicity', '')
                am_str = f"{am:.2f}" if isinstance(am, (int, float)) and am else str(am) if am else '—'
                ev_data.append([
                    self._safe(str(ev.get('gene', '—'))),
                    self._safe(str(ev.get('rsid', '—'))),
                    self._safe(str(ev.get('consequence', '—'))),
                    self._safe(str(ev.get('impact', '—'))),
                    score_str,
                    am_str,
                ])

            ev_table = Table(ev_data, colWidths=[0.8*inch, 1.2*inch, 1.2*inch, 0.8*inch, 0.7*inch, 1.1*inch])
            ev_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), Colors.ACCENT_BLUE),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f0f7ff')]),
                ('BOX', (0, 0), (-1, -1), 1, Colors.ACCENT_BLUE),
                ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#cce0ff')),
                ('TOPPADDING', (0, 0), (-1, -1), 5),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
            ]))
            elements.append(KeepTogether([ev_header, ev_table]))
            elements.append(Spacer(1, 10))

        # Knowledge connections
        connections = self.context.get('knowledge_connections', [])
        if connections:
            kc_header = Paragraph(f"Knowledge Connections ({len(connections)} genes)", self.styles['SubSection'])

            kc_data = [['Gene', 'Protein', 'Pathway', 'Drugs', 'Druggable']]
            for kc in connections[:8]:
                drugs = kc.get('drugs', [])
                drug_str = ', '.join(drugs[:2]) + ('...' if len(drugs) > 2 else '') if drugs else '—'
                kc_data.append([
                    self._safe(str(kc.get('gene', ''))),
                    self._safe(str(kc.get('protein', ''))[:30]),
                    self._safe(str(kc.get('pathway', ''))[:25]),
                    self._safe(drug_str[:30]),
                    'YES' if kc.get('druggable') else 'NO',
                ])

            kc_table = Table(kc_data, colWidths=[0.7*inch, 1.5*inch, 1.3*inch, 1.5*inch, 0.8*inch])
            kc_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), Colors.ACCENT_PURPLE),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f8fc')]),
                ('BOX', (0, 0), (-1, -1), 1, Colors.ACCENT_PURPLE),
                ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e0e0f0')),
                ('TOPPADDING', (0, 0), (-1, -1), 5),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
            ]))
            elements.append(KeepTogether([kc_header, kc_table]))

        return elements

    @staticmethod
    def _safe(text: str) -> str:
        """Escape XML special characters for ReportLab paragraphs."""
        return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

    @staticmethod
    def _num(value, places: int) -> str:
        """Format a number, or say it was not computed. Never substitutes a default."""
        if isinstance(value, bool) or value is None:
            return 'not computed'
        if isinstance(value, (int, float)):
            return f"{value:.{places}f}"
        return str(value)

    @staticmethod
    def _dock_cell(candidate: dict) -> str:
        """Render a docking result in its own units, or say there was none.

        The old table printed every value as kcal/mol. DiffDock reports a pose confidence, not
        an affinity, so a confidence rendered as kcal/mol would be a unit the method never had.
        """
        score_type = candidate.get('score_type')
        if score_type == 'diffdock_confidence':
            conf = candidate.get('best_pose_confidence', candidate.get('pose_confidence'))
            return f"pose conf. {conf:.3f}" if isinstance(conf, (int, float)) else 'not docked'
        score = candidate.get('docking_score')
        if isinstance(score, (int, float)):
            return f"{score:.2f} kcal/mol"
        conf = candidate.get('best_pose_confidence', candidate.get('pose_confidence'))
        if isinstance(conf, (int, float)) and conf:
            return f"pose conf. {conf:.3f}"
        return 'not docked'

    @staticmethod
    def _descriptors(smiles: str):
        """RDKit descriptors for a SMILES, or None when they cannot be computed."""
        if not smiles:
            return None
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
        except ImportError:
            return None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return {
            'mw': round(Descriptors.MolWt(mol), 2),
            'logp': round(Descriptors.MolLogP(mol), 2),
            'hbd': Lipinski.NumHDonors(mol),
            'hba': Lipinski.NumHAcceptors(mol),
        }

    def _create_structure_section(self):
        """Create structural biology section with Cryo-EM images."""
        elements = []

        elements.append(PageBreak())
        elements.append(Paragraph("3. STRUCTURAL EVIDENCE", self.styles['SectionHeader']))

        elements.append(Paragraph(
            "High-resolution Cryo-EM and X-ray structures of VCP/p97 provide molecular "
            "templates for structure-based drug design. These structures reveal the "
            "validated D2 ATPase binding pocket, targeted by ATP-competitive inhibitors such as CB-5083.",
            self.styles['BodyJustified']
        ))
        elements.append(Spacer(1, 15))

        # Structure gallery with images
        elements.append(Paragraph("VCP Structure Gallery", self.styles['SubSection']))

        # Create 2x2 grid of structure images
        structure_images = []
        for struct in self.structures[:4]:
            pdb_id = struct['structure_id'].replace('PDB:', '')
            img_path = self._get_structure_image(pdb_id)

            if img_path and img_path.exists():
                try:
                    img = Image(str(img_path), width=2.8*inch, height=2.1*inch)
                    caption = Paragraph(
                        f"<b>{pdb_id}</b> | {struct['method']} | {struct['resolution']}<br/>"
                        f"<font size='8'>{struct['conformation']}</font>",
                        self.styles['Caption']
                    )
                    structure_images.append([img, caption])
                except Exception as e:
                    print(f"Could not load image for {pdb_id}: {e}")
                    structure_images.append([
                        Paragraph(f"<b>{pdb_id}</b><br/>{struct['method']}<br/>{struct['resolution']}",
                                  self.styles['BodyJustified']),
                        Paragraph(struct['conformation'], self.styles['Caption'])
                    ])
            else:
                structure_images.append([
                    Paragraph(f"<b>{pdb_id}</b><br/>{struct['method']}<br/>{struct['resolution']}",
                              self.styles['BodyJustified']),
                    Paragraph(struct['conformation'], self.styles['Caption'])
                ])

        # Arrange in 2x2 grid
        if len(structure_images) >= 4:
            grid_data = [
                [structure_images[0], structure_images[1]],
                [structure_images[2], structure_images[3]],
            ]

            # Flatten the nested structure
            flat_grid = []
            for row in grid_data:
                flat_row = []
                for cell in row:
                    cell_content = []
                    for item in cell:
                        cell_content.append(item)
                    flat_row.append(cell_content)
                flat_grid.append(flat_row)

            for row in flat_grid:
                row_table_data = [[row[0][0], row[1][0]], [row[0][1], row[1][1]]]
                row_table = Table(row_table_data, colWidths=[3.1*inch, 3.1*inch])
                row_table.setStyle(TableStyle([
                    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                    ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                    ('BOX', (0, 0), (-1, -1), 1, colors.HexColor('#e0e0e0')),
                    ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#fafafa')),
                    ('TOPPADDING', (0, 0), (-1, -1), 8),
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
                ]))
                elements.append(row_table)
                elements.append(Spacer(1, 10))

        # Primary docking structure.
        #
        # This block used to assert a "2.3 A crystal structure of VCP bound to CB-5083" for PDB
        # 5FTK (the orchestrator was corrected to call the VCP entries cryo-EM on 2026-11-03, so
        # the two files contradicted each other) followed by a binding-site table: key residues
        # ALA464/GLY479/ASP320/GLY215, pocket volume ~450 A^3, druggability score 0.92. No pocket
        # detection, volume calculation or druggability model runs anywhere in this repo -- that
        # is the same class of invented binding-site evidence that was deleted from
        # modules/drug_discovery.nf, and two of the four residues are N/D1-domain numbering
        # presented as D2 pocket residues. Removed outright.
        elements.append(Spacer(1, 10))
        meta = self._pipeline_meta()
        template = meta.get("pdb_id")
        docking_header = Paragraph(
            f"Primary Docking Template: PDB {self._safe(str(template))}" if template
            else "Primary Docking Template",
            self.styles['SubSection'],
        )
        docking_desc = Paragraph(
            (f"PDB {self._safe(str(template))} was the receptor supplied to docking for this run, "
             "as deposited in the RCSB PDB. "
             if template else
             "No docking template was supplied to this report. ")
            + "Its experimental method and resolution are the deposited values shown in the "
              "gallery above; this report performs no pocket detection and reports no pocket "
              "volume, key-residue list or druggability score, because it computes none.",
            self.styles['BodyJustified']
        )
        elements.append(KeepTogether([docking_header, docking_desc]))

        return elements

    def _create_molecules_section(self):
        """Create drug candidate section from the supplied run, or refuse.

        This method was the heart of the fabrication. It built "Top Ranked Drug Candidates" from
        four literal dicts -- VCP-AI-001..004 with `docking` -8.62 to -10.95 kcal/mol, `qed`,
        `mw`, `logp` and `score` -- under the prose "NVIDIA BioNeMo MolMIM generated novel
        molecules ... docked against VCP using DiffDock", and closed with "All candidates satisfy
        Lipinski's Rule of Five and show favorable ADMET predictions". The file calls no NIM, no
        docking code, no RDKit descriptor and no ADMET predictor. The table was not even
        self-consistent with the repo's own composite formula: candidate #4 had the best docking
        value and the worst composite score. It now renders the run it is given, and raises when
        it is given none.
        """
        self._require_results("Section 4 (generated drug candidates)")

        elements = []
        meta = self._pipeline_meta()
        candidates = self._candidates()

        elements.append(PageBreak())
        elements.append(Paragraph("4. GENERATED DRUG CANDIDATES", self.styles['SectionHeader']))

        # Name only the tools the run recorded. An unrecorded step is reported as unrecorded.
        generator_name = meta.get("generator") or meta.get("nim_mode") or "not recorded"
        docking_method = meta.get("docking_method")
        n_generated = meta.get("num_generated", len(candidates))
        n_passed = meta.get("num_passed_qc")
        n_docked = meta.get("num_docked", 0)

        prose = (
            f"Candidates below come from the run supplied to this report. "
            f"Molecule generation: <b>{self._safe(str(generator_name))}</b>; "
            f"{self._safe(str(n_generated))} molecule(s) generated"
            + (f", {self._safe(str(n_passed))} passed chemistry QC" if n_passed is not None else "")
            + ". "
        )
        if n_docked and docking_method:
            prose += (
                f"Docking: <b>{self._safe(str(docking_method))}</b> against "
                f"{self._safe(str(meta.get('pdb_id', 'the supplied receptor')))}, "
                f"{self._safe(str(n_docked))} molecule(s) docked. A pose confidence is a "
                "probability that the pose is near-native; it is not a binding affinity, so no "
                "potency is implied."
            )
        else:
            prose += (
                "<b>No docking was performed in this run</b>, so no binding score of any kind is "
                "reported below and the ranking carries no binding term."
            )
        elements.append(Paragraph(prose, self.styles['BodyJustified']))
        elements.append(Spacer(1, 10))

        # Reference compound with image — identity and properties come from the run.
        seed_name = meta.get("seed_compound") or "not recorded"
        seed_smiles = meta.get("seed_smiles") or ""
        elements.append(Paragraph(
            f"Seed Compound: {self._safe(str(seed_name))}", self.styles['SubSection']))

        mol_img_path = (self._download_molecule_image(seed_smiles, "seed")
                        if seed_smiles else None)

        ref_content = []
        if mol_img_path and mol_img_path.exists():
            try:
                mol_img = Image(str(mol_img_path), width=2*inch, height=2*inch)
                ref_content.append([mol_img])
            except Exception:
                pass

        # The seed's MW/LogP/HBD/HBA used to be hard-coded (484.6 / 4.92 / 2 / 6). Compute them
        # from the SMILES, or report them as not computed; never print a remembered number.
        ref_props = [['Property', 'Value']]
        seed_desc = self._descriptors(seed_smiles)
        if seed_desc:
            ref_props += [
                ['MW', f"{seed_desc['mw']} Da"],
                ['LogP', f"{seed_desc['logp']}"],
                ['HBD', f"{seed_desc['hbd']}"],
                ['HBA', f"{seed_desc['hba']}"],
            ]
        else:
            ref_props.append(['Descriptors', 'not computed (no SMILES or RDKit unavailable)'])

        ref_table = Table(ref_props, colWidths=[1.5*inch, 2*inch])
        ref_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), Colors.NVIDIA_GREEN),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8fff0')]),
            ('BOX', (0, 0), (-1, -1), 1, Colors.NVIDIA_GREEN),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#c8e8a0')),
            ('TOPPADDING', (0, 0), (-1, -1), 5),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
        ]))

        if ref_content:
            layout_data = [[ref_content[0][0], ref_table]]
            layout = Table(layout_data, colWidths=[2.5*inch, 3.5*inch])
            layout.setStyle(TableStyle([
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ]))
            elements.append(layout)
        else:
            elements.append(ref_table)

        if seed_smiles:
            elements.append(Spacer(1, 8))
            elements.append(Paragraph(
                f"<font face='Courier' size='7'>SMILES: {self._safe(seed_smiles)}</font>",
                self.styles['SmallText']
            ))
        elements.append(Spacer(1, 15))

        # Top candidates — every cell from the supplied run.
        elements.append(Paragraph("Top Ranked Drug Candidates", self.styles['SubSection']))

        cand_header = ['Rank', 'ID', 'Docking', 'QED', 'MW\n(Da)', 'LogP', 'Composite', 'Provenance']
        cand_data = [cand_header]

        for i, c in enumerate(candidates[:10], start=1):
            props = c.get('properties') or {}
            cand_data.append([
                f"#{c.get('rank', i)}",
                self._safe(str(c.get('id') or c.get('molecule_id') or c.get('name') or f'mol-{i}')),
                self._safe(self._dock_cell(c)),
                self._safe(self._num(props.get('qed', c.get('qed_score')), 3)),
                self._safe(self._num(props.get('mw', props.get('molecular_weight')), 1)),
                self._safe(self._num(props.get('logp', props.get('logP')), 2)),
                self._safe(self._num(c.get('composite_score', c.get('score')), 4)),
                self._safe(str(props.get('provenance') or c.get('provenance') or 'not recorded')),
            ])

        cand_table = Table(cand_data, colWidths=[0.5*inch, 1.0*inch, 1.0*inch, 0.55*inch,
                                                 0.6*inch, 0.55*inch, 0.75*inch, 1.25*inch])
        cand_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), Colors.NVIDIA_GREEN),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            # Highlight rank 1
            ('BACKGROUND', (0, 1), (-1, 1), colors.HexColor('#e8f5e0')),
            ('FONTNAME', (0, 1), (-1, 1), 'Helvetica-Bold'),
            ('ROWBACKGROUNDS', (0, 2), (-1, -1), [colors.white, colors.HexColor('#f8f8f8')]),
            ('BOX', (0, 0), (-1, -1), 2, Colors.NVIDIA_GREEN),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#c8e8a0')),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ]))
        elements.append(cand_table)
        elements.append(Spacer(1, 15))

        # Drug-likeness assessment.
        #
        # This block used to assert "All candidates satisfy Lipinski's Rule of Five and show
        # favorable ADMET predictions" and then tabulate four PASS rows for VCP-AI-001 from the
        # same hard-coded numbers. No ADMET predictor is invoked in this file at all (the engine
        # has one, src/admet_predictor.py, and this report never calls it), so the ADMET clause
        # is gone, and the Lipinski table is built from the top candidate's own descriptors.
        top = candidates[0]
        top_props = top.get('properties') or {}
        top_id = str(top.get('id') or top.get('molecule_id') or top.get('name') or 'rank 1')
        lip_header = Paragraph("Drug-Likeness Assessment", self.styles['SubSection'])
        lip_desc = Paragraph(
            "Lipinski Rule-of-Five terms for the top-ranked candidate, from the descriptors "
            "supplied with the run. No ADMET model was run for this report, so no ADMET "
            "prediction is reported.",
            self.styles['BodyJustified']
        )

        def _rule(label, threshold, value, ok):
            if value is None:
                return [label, threshold, 'not computed', 'NOT ASSESSED']
            return [label, threshold, str(value), 'PASS' if ok(value) else 'FAIL']

        mw_v = top_props.get('mw', top_props.get('molecular_weight'))
        logp_v = top_props.get('logp', top_props.get('logP'))
        lipinski_data = [
            ['Rule', 'Threshold', self._safe(top_id), 'Status'],
            _rule('Molecular Weight', '<= 500 Da', mw_v, lambda v: v <= 500),
            _rule('LogP', '<= 5', logp_v, lambda v: v <= 5),
            _rule('H-Bond Donors', '<= 5', top_props.get('hbd'), lambda v: v <= 5),
            _rule('H-Bond Acceptors', '<= 10', top_props.get('hba'), lambda v: v <= 10),
        ]

        lip_table = Table(lipinski_data, colWidths=[1.5*inch, 1.3*inch, 1.3*inch, 1.2*inch])
        lip_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), Colors.ACCENT_CYAN),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
            ('TEXTCOLOR', (3, 1), (3, -1), Colors.STATUS_SUCCESS),
            ('FONTNAME', (3, 1), (3, -1), 'Helvetica-Bold'),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f0ffff')]),
            ('BOX', (0, 0), (-1, -1), 1, Colors.ACCENT_CYAN),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#b0e0e0')),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
        ]))
        # Keep drug-likeness section together
        elements.append(KeepTogether([lip_header, lip_desc, Spacer(1, 8), lip_table]))

        return elements

    def _create_summary_section(self):
        """Create executive summary from the supplied run, or refuse.

        The four "key achievements" were literals: an AlphaMissense score of 0.89 for a variant
        nothing called, "VCP/p97 confirmed", "4 high-resolution structures (2.3-3.2 A)" and
        "4 drug candidates with docking scores -8.26 to -10.95 kcal/mol" -- the fabricated table
        from section 4, quoted back as a result. The "Pipeline Performance" table that followed
        gave stage timings (Genomics 120-240 min on Parabricks 4.6, Molecule Generation 2-5 min
        on BioNeMo MolMIM, Docking 5-10 min on DiffDock, total < 5 hours) that nothing timed and
        that the sibling generator contradicted (30-110 min for the same Parabricks stage).
        """
        self._require_results("Section 5 (executive summary)")

        elements = []
        meta = self._pipeline_meta()
        candidates = self._candidates()

        elements.append(PageBreak())
        elements.append(Paragraph("5. RUN SUMMARY", self.styles['SectionHeader']))

        # What this run actually did, counted from the run itself.
        facts = [
            ("Molecules generated", str(meta.get("num_generated", len(candidates)))),
            ("Passed chemistry QC", str(meta.get("num_passed_qc", "not recorded"))),
            ("Molecules docked", str(meta.get("num_docked", 0))),
            ("Generator", str(meta.get("generator") or meta.get("nim_mode") or "not recorded")),
            ("Docking method", str(meta.get("docking_method") or "no docking performed")),
        ]
        for title, text in facts:
            elements.append(Paragraph(
                f"<font color='#76B900'>●</font> <b>{self._safe(title)}:</b> {self._safe(text)}",
                self.styles['BodyJustified']
            ))

        # Only report timings the run measured. Nothing is filled in from a remembered range.
        timings = {k: v for k, v in (
            ("Molecule generation", meta.get("generation_time_sec")),
            ("Docking", meta.get("docking_time_sec")),
            ("Total", meta.get("total_time_sec")),
        ) if v is not None}
        if timings:
            elements.append(Spacer(1, 10))
            elements.append(Paragraph("Measured wall-clock (this run)", self.styles['SubSection']))
            for label, secs in timings.items():
                elements.append(Paragraph(
                    f"{self._safe(label)}: {self._safe(str(secs))} s",
                    self.styles['BodyJustified']
                ))
        else:
            elements.append(Spacer(1, 10))
            elements.append(Paragraph(
                "<b>No stage timings were recorded for this run</b>, so none are reported.",
                self.styles['BodyJustified']
            ))

        elements.append(Spacer(1, 20))
        return elements

    def _create_next_steps_section(self):
        """Laboratory follow-up suggestions. Not a result — labelled as such."""
        elements = []
        elements.append(Paragraph("Suggested Laboratory Follow-up", self.styles['SubSection']))
        elements.append(Paragraph(
            "Generic wet-lab follow-up for any VCP-directed series. These are suggestions, not "
            "findings, and are not derived from this run's numbers.",
            self.styles['SmallText']
        ))

        steps = [
            "Synthesize the top candidates for biochemical validation",
            "Evaluate VCP ATPase inhibition in enzymatic assays",
            "Assess blood-brain barrier permeability for CNS penetration",
            "Profile selectivity against related AAA+ ATPases",
            "Evaluate in cellular models of VCP-associated disease",
        ]
        for i, step in enumerate(steps, 1):
            elements.append(Paragraph(
                f"<font color='#76B900'><b>{i}.</b></font> {step}",
                self.styles['BodyJustified']
            ))
        return elements

    def _create_footer(self):
        """Create document footer."""
        elements = []

        elements.append(Spacer(1, 30))
        elements.append(GradientRect(6.5*inch, 2, Colors.NVIDIA_GREEN, Colors.NVIDIA_GREEN_LIGHT))
        elements.append(Spacer(1, 15))

        # The platform line used to credit Parabricks, BioNeMo NIM, Milvus and Claude on every
        # report, whatever had actually run. Credit only what the run recorded.
        meta = self._pipeline_meta()
        ran = [str(v) for v in (meta.get("generator") or meta.get("nim_mode"),
                                meta.get("docking_method")) if v]
        tools = " | ".join(ran) if ran else "no tool run recorded"

        footer = Paragraph(
            f"<b>HCLS AI Factory</b> | Precision Medicine to Drug Discovery<br/>"
            f"Generated: {datetime.now().strftime('%B %d, %Y at %H:%M')}<br/>"
            f"Tools recorded for this run: {self._safe(tools)}<br/>"
            f"<br/>"
            f"Research use only — not for clinical decision-making.",
            self.styles['SmallText']
        )
        footer.style.alignment = TA_CENTER
        elements.append(footer)

        return elements

    def generate(self):
        """Generate the complete PDF report."""
        doc = SimpleDocTemplate(
            str(self.output_path),
            pagesize=letter,
            rightMargin=0.6*inch,
            leftMargin=0.6*inch,
            topMargin=0.6*inch,
            bottomMargin=0.6*inch
        )

        elements = []

        # Build report sections
        elements.extend(self._create_header())
        elements.extend(self._create_pipeline_flow())
        elements.extend(self._create_variant_section())
        elements.extend(self._create_target_section())
        if self.context:
            elements.extend(self._create_rag_query_section())
        elements.extend(self._create_structure_section())
        elements.extend(self._create_molecules_section())
        elements.extend(self._create_summary_section())
        elements.extend(self._create_next_steps_section())
        elements.extend(self._create_footer())

        # Build PDF
        doc.build(elements)
        print(f"Report generated: {self.output_path}")
        return self.output_path


def main():
    """Generate enhanced VCP Drug Candidate Report."""
    import argparse
    import shutil

    parser = argparse.ArgumentParser(description='Generate VCP Drug Candidate PDF Report')
    parser.add_argument('--context', type=str, help='Path to report_context.json from RAG chat')
    parser.add_argument('--results', type=str, required=True,
                        help='Path to a pipeline results JSON (e.g. outputs/cloud_nim_results.json '
                             'written by run_cloud_nim_report.py). Required: this generator '
                             'reports a run, it does not invent one.')
    args = parser.parse_args()

    # Load context if provided
    context = None
    if args.context and Path(args.context).exists():
        with open(args.context) as f:
            context = json.load(f)

    results_path = Path(args.results)
    if not results_path.exists():
        raise SystemExit(
            f"--results {args.results} does not exist. Run the pipeline first "
            "(python run_cloud_nim_report.py) and point --results at its JSON output. "
            "No report is written without a run to report."
        )
    with open(results_path) as f:
        results = json.load(f)

    generator = VCPReportGeneratorEnhanced(
        output_path="outputs/VCP_Drug_Candidate_Report.pdf",
        context=context,
        results=results,
    )
    output_file = generator.generate()

    # Copy to landing page
    landing_page_path = Path("../landing-page/static/VCP_Drug_Candidate_Report.pdf")
    if landing_page_path.parent.exists():
        shutil.copy(output_file, landing_page_path)
        print(f"Copied to: {landing_page_path}")

    # Write metadata for landing page freshness check
    meta_path = Path("../landing-page/static/report_meta.json")
    if meta_path.parent.exists():
        meta = {
            'generated_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'has_context': context is not None,
            'query': context.get('query', '')[:100] if context else '',
        }
        with open(meta_path, 'w') as f:
            json.dump(meta, f)

    print(f"Report generated: {output_file}")


if __name__ == "__main__":
    main()
