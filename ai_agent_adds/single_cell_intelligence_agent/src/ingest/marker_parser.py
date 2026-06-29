"""Marker gene parser for the Single-Cell Intelligence Agent.

Seeds 75 marker gene records from CellMarker and PanglaoDB databases,
covering canonical markers for immune, epithelial, stromal, neural, and
stem cell populations used in single-cell RNA-seq cell type annotation.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: 75 MARKER GENE RECORDS FROM CellMarker / PanglaoDB
# ===================================================================

MARKER_GENE_RECORDS: List[Dict[str, Any]] = [
    # --- T cell markers ---
    {"gene": "CD3D", "cell_types": ["T cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "CD3 delta chain, part of the TCR-CD3 complex. Pan-T cell marker essential for T cell development and signaling."},
    {"gene": "CD3E", "cell_types": ["T cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "CD3 epsilon chain. Required for surface expression of TCR complex. Used as a pan-T cell marker."},
    {"gene": "CD4", "cell_types": ["CD4+ T cell", "Monocyte"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "Glycoprotein functioning as a co-receptor for MHC class II. Defines helper T cell lineage."},
    {"gene": "CD8A", "cell_types": ["CD8+ T cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "CD8 alpha chain, co-receptor for MHC class I. Defines cytotoxic T cell lineage."},
    {"gene": "FOXP3", "cell_types": ["Regulatory T cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Forkhead box P3 transcription factor, master regulator of regulatory T cell development and function."},
    {"gene": "GZMB", "cell_types": ["CD8+ T cell", "NK cell"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "Granzyme B serine protease. Cytotoxic effector molecule in CD8+ T cells and NK cells."},
    {"gene": "PRF1", "cell_types": ["CD8+ T cell", "NK cell"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "Perforin 1 pore-forming protein. Creates pores in target cell membranes for granzyme delivery."},
    {"gene": "IL7R", "cell_types": ["T cell", "ILC"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "IL-7 receptor alpha chain (CD127). Expressed on naive and memory T cells. Low on Tregs."},
    # --- B cell markers ---
    {"gene": "CD19", "cell_types": ["B cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "B-lymphocyte antigen CD19. Pan-B cell marker from pro-B to memory B cell stages. Target of CAR-T therapy."},
    {"gene": "MS4A1", "cell_types": ["B cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Membrane-spanning 4A1, also known as CD20. Target of rituximab. Lost upon plasma cell differentiation."},
    {"gene": "CD79A", "cell_types": ["B cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "CD79a immunoglobulin-associated alpha. Component of the B cell receptor signaling complex."},
    {"gene": "JCHAIN", "cell_types": ["Plasma cell"], "specificity": "high", "source_db": "panglaodb",
     "description": "Immunoglobulin J chain. Essential for IgM pentamer and IgA dimer assembly in plasma cells."},
    {"gene": "SDC1", "cell_types": ["Plasma cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Syndecan-1 (CD138). Cell surface proteoglycan highly expressed on plasma cells. Clinical flow cytometry marker."},
    # --- NK cell markers ---
    {"gene": "NCAM1", "cell_types": ["NK cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Neural cell adhesion molecule 1 (CD56). Principal surface marker for NK cell identification."},
    {"gene": "KLRD1", "cell_types": ["NK cell"], "specificity": "high", "source_db": "panglaodb",
     "description": "Killer cell lectin-like receptor D1 (CD94). Forms heterodimers with NKG2 family members on NK cells."},
    {"gene": "NKG7", "cell_types": ["NK cell", "CD8+ T cell"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "Natural killer cell granule protein 7. Cytotoxic granule membrane protein in NK and cytotoxic T cells."},
    # --- Myeloid markers ---
    {"gene": "CD14", "cell_types": ["Monocyte", "Macrophage"], "specificity": "high", "source_db": "cellmarker",
     "description": "Monocyte differentiation antigen CD14. LPS co-receptor. Classical monocyte marker (CD14++CD16-)."},
    {"gene": "CD68", "cell_types": ["Macrophage"], "specificity": "high", "source_db": "cellmarker",
     "description": "Macrosialin/CD68 glycoprotein. Widely used macrophage marker expressed in lysosomes and cell surface."},
    {"gene": "CD163", "cell_types": ["Macrophage"], "specificity": "high", "source_db": "cellmarker",
     "description": "Hemoglobin scavenger receptor. M2 macrophage marker associated with anti-inflammatory polarization."},
    {"gene": "ITGAX", "cell_types": ["Dendritic cell", "Macrophage"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "Integrin alpha X (CD11c). Key dendritic cell marker. Also expressed on macrophages and some monocytes."},
    {"gene": "CLEC9A", "cell_types": ["cDC1"], "specificity": "high", "source_db": "cellmarker",
     "description": "C-type lectin domain 9A. Specific marker for conventional type 1 dendritic cells (cDC1) that cross-present antigens."},
    {"gene": "S100A8", "cell_types": ["Monocyte", "Neutrophil"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "S100 calcium-binding protein A8 (calgranulin A). Inflammatory marker in monocytes and neutrophils."},
    {"gene": "FCGR3A", "cell_types": ["NK cell", "Monocyte"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "Fc gamma receptor IIIa (CD16). Expressed on non-classical monocytes and CD56dim NK cells."},
    {"gene": "CSF3R", "cell_types": ["Neutrophil"], "specificity": "high", "source_db": "cellmarker",
     "description": "Colony stimulating factor 3 receptor (G-CSFR). Key neutrophil lineage marker and G-CSF signaling receptor."},
    {"gene": "KIT", "cell_types": ["Mast cell", "HSC"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "KIT proto-oncogene (CD117). Stem cell factor receptor. Marker for mast cells and hematopoietic stem cells."},
    {"gene": "TPSAB1", "cell_types": ["Mast cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Tryptase alpha/beta 1. Mast cell-specific serine protease. Diagnostic marker for systemic mastocytosis."},
    # --- Epithelial markers ---
    {"gene": "EPCAM", "cell_types": ["Epithelial cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Epithelial cell adhesion molecule. Pan-epithelial marker used in flow cytometry and CTC detection."},
    {"gene": "KRT18", "cell_types": ["Epithelial cell"], "specificity": "high", "source_db": "panglaodb",
     "description": "Keratin 18 intermediate filament. Simple epithelium marker. Expressed in most adenocarcinomas."},
    {"gene": "CDH1", "cell_types": ["Epithelial cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "E-cadherin. Calcium-dependent cell adhesion. Loss marks epithelial-mesenchymal transition (EMT)."},
    {"gene": "SFTPC", "cell_types": ["AT2 cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Surfactant protein C. Specific marker for alveolar type II cells. Mutations cause interstitial lung disease."},
    {"gene": "MUC2", "cell_types": ["Goblet cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Mucin 2. Gel-forming mucin produced by intestinal goblet cells. Forms the protective mucus barrier."},
    # --- Stromal/Mesenchymal markers ---
    {"gene": "COL1A1", "cell_types": ["Fibroblast", "Osteoblast"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "Collagen type I alpha 1. Major structural protein in fibroblasts. Mutations cause osteogenesis imperfecta."},
    {"gene": "DCN", "cell_types": ["Fibroblast"], "specificity": "high", "source_db": "panglaodb",
     "description": "Decorin proteoglycan. Fibroblast marker that binds TGF-beta and inhibits fibrosis."},
    {"gene": "FAP", "cell_types": ["Cancer-associated fibroblast"], "specificity": "high", "source_db": "cellmarker",
     "description": "Fibroblast activation protein. Serine protease enriched in CAFs and reactive stroma. Therapeutic target."},
    {"gene": "PDGFRA", "cell_types": ["Fibroblast", "MSC"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "PDGF receptor alpha. Marks fibroblast and mesenchymal progenitor populations. Target in GIST therapy."},
    {"gene": "PDGFRB", "cell_types": ["Pericyte", "Fibroblast"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "PDGF receptor beta. Pericyte marker involved in vascular stability and recruitment."},
    {"gene": "ACTA2", "cell_types": ["Smooth muscle cell", "Myofibroblast", "Pericyte"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "Alpha smooth muscle actin. Marks smooth muscle cells and activated myofibroblasts."},
    # --- Endothelial markers ---
    {"gene": "PECAM1", "cell_types": ["Endothelial cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Platelet endothelial cell adhesion molecule 1 (CD31). Pan-endothelial marker. Mediates leukocyte transmigration."},
    {"gene": "VWF", "cell_types": ["Endothelial cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Von Willebrand factor. Endothelial-specific glycoprotein involved in hemostasis. Stored in Weibel-Palade bodies."},
    {"gene": "CDH5", "cell_types": ["Endothelial cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "VE-cadherin. Endothelial-specific adherens junction protein. Critical for vascular permeability control."},
    # --- Stem cell markers ---
    {"gene": "CD34", "cell_types": ["HSC", "Endothelial progenitor"], "specificity": "high", "source_db": "cellmarker",
     "description": "Hematopoietic progenitor cell antigen CD34. Marks HSCs and endothelial progenitors. Used for transplant selection."},
    {"gene": "POU5F1", "cell_types": ["Embryonic stem cell", "iPSC"], "specificity": "high", "source_db": "cellmarker",
     "description": "OCT4 transcription factor. Core pluripotency factor in ESCs and iPSCs. Silenced upon differentiation."},
    {"gene": "NANOG", "cell_types": ["Embryonic stem cell", "iPSC"], "specificity": "high", "source_db": "cellmarker",
     "description": "Homeobox protein NANOG. Maintains pluripotency alongside OCT4 and SOX2. Lost during differentiation."},
    {"gene": "SOX2", "cell_types": ["Stem cell", "Neural progenitor"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "SRY-box 2 transcription factor. Pluripotency factor and neural progenitor marker."},
    {"gene": "LGR5", "cell_types": ["Intestinal stem cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Leucine-rich repeat G-protein coupled receptor 5. Wnt target gene marking intestinal crypt stem cells."},
    # --- Neural markers ---
    {"gene": "RBFOX3", "cell_types": ["Neuron"], "specificity": "high", "source_db": "cellmarker",
     "description": "RNA-binding Fox protein 3 (NeuN). Post-mitotic neuron-specific nuclear marker used in immunohistochemistry."},
    {"gene": "GFAP", "cell_types": ["Astrocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "Glial fibrillary acidic protein. Intermediate filament in astrocytes. Upregulated in reactive astrogliosis."},
    {"gene": "MBP", "cell_types": ["Oligodendrocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "Myelin basic protein. Major component of CNS myelin. Target of autoimmune attack in multiple sclerosis."},
    {"gene": "TMEM119", "cell_types": ["Microglia"], "specificity": "high", "source_db": "cellmarker",
     "description": "Transmembrane protein 119. Microglia-specific marker distinguishing brain-resident microglia from infiltrating macrophages."},
    {"gene": "OLIG2", "cell_types": ["Oligodendrocyte", "OPC"], "specificity": "high", "source_db": "cellmarker",
     "description": "Oligodendrocyte transcription factor 2. Essential for oligodendrocyte lineage specification."},
    # --- Other tissue markers ---
    {"gene": "ALB", "cell_types": ["Hepatocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "Albumin. Major plasma protein produced exclusively by hepatocytes. Gold standard liver cell marker."},
    {"gene": "ADIPOQ", "cell_types": ["Adipocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "Adiponectin. Adipocyte-specific hormone with insulin-sensitizing and anti-inflammatory properties."},
    {"gene": "MITF", "cell_types": ["Melanocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "Melanocyte inducing transcription factor. Master regulator of melanocyte development and melanoma."},
    {"gene": "RUNX2", "cell_types": ["Osteoblast"], "specificity": "high", "source_db": "cellmarker",
     "description": "Runt-related transcription factor 2. Master regulator of osteoblast differentiation and bone formation."},
    {"gene": "SOX9", "cell_types": ["Chondrocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "SRY-box 9 transcription factor. Master regulator of chondrocyte differentiation and cartilage formation."},
    # --- Gamma-delta T / unconventional T cell markers ---
    {"gene": "TRGV9", "cell_types": ["Gamma-delta T cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "T cell receptor gamma variable 9. Defines the Vgamma9 chain of the predominant Vgamma9Vdelta2 T cell subset in peripheral blood."},
    {"gene": "SLC4A10", "cell_types": ["MAIT cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Solute carrier family 4 member 10 (sodium bicarbonate transporter). Highly specific surface marker for MAIT cells in scRNA-seq."},
    # --- Dendritic cell subtype markers ---
    {"gene": "CLEC4C", "cell_types": ["Plasmacytoid dendritic cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "C-type lectin domain family 4 member C (BDCA-2/CD303). Specific marker for plasmacytoid dendritic cells."},
    {"gene": "XCR1", "cell_types": ["cDC1"], "specificity": "high", "source_db": "cellmarker",
     "description": "X-C motif chemokine receptor 1. Specific marker for cDC1 cells. Ligand XCL1 produced by CD8+ T and NK cells recruits cDC1 to tumors."},
    # --- ILC markers ---
    {"gene": "RORC", "cell_types": ["ILC3", "Th17"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "RAR-related orphan receptor C (RORgamma-t). Master transcription factor for ILC3 and Th17 cell differentiation."},
    # --- Megakaryocyte / platelet markers ---
    {"gene": "ITGA2B", "cell_types": ["Megakaryocyte", "Platelet"], "specificity": "high", "source_db": "cellmarker",
     "description": "Integrin alpha-2b (CD41/GPIIb). Surface glycoprotein marking megakaryocyte lineage. Forms GPIIb/IIIa complex essential for platelet aggregation."},
    # --- Erythroid markers ---
    {"gene": "GYPA", "cell_types": ["Erythroid progenitor", "Erythrocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "Glycophorin A (CD235a). Major sialoglycoprotein on erythrocyte surface. Definitive erythroid lineage marker."},
    # --- Pericyte / vascular markers ---
    {"gene": "RGS5", "cell_types": ["Pericyte"], "specificity": "high", "source_db": "panglaodb",
     "description": "Regulator of G-protein signaling 5. Highly specific pericyte marker involved in vascular remodeling and angiogenesis."},
    # --- Schwann cell / PNS markers ---
    {"gene": "MPZ", "cell_types": ["Schwann cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Myelin protein zero. Major structural protein of peripheral nerve myelin produced by Schwann cells. Mutations cause Charcot-Marie-Tooth disease."},
    # --- Podocyte / kidney markers ---
    {"gene": "NPHS1", "cell_types": ["Podocyte"], "specificity": "high", "source_db": "cellmarker",
     "description": "Nephrin. Key component of the slit diaphragm in kidney podocytes. Mutations cause congenital nephrotic syndrome."},
    # --- Goblet cell / secretory markers ---
    {"gene": "MUC2", "cell_types": ["Goblet cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Mucin 2. Gel-forming mucin produced by intestinal goblet cells. Forms the protective mucus barrier in the colon."},
    # --- Stromal / fibroblast markers ---
    {"gene": "POSTN", "cell_types": ["Cancer-associated fibroblast", "Fibroblast"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "Periostin extracellular matrix protein. Marks activated fibroblasts and CAFs. Promotes tumor invasion and metastatic niche formation."},
    # --- T cell state markers ---
    {"gene": "TCF7", "cell_types": ["Naive T cell", "Progenitor exhausted T cell"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "T cell factor 7 (TCF-1). Marks naive and stem-like/progenitor exhausted T cells. TCF7+ CD8 T cells sustain anti-tumor responses during checkpoint blockade."},
    {"gene": "CXCR6", "cell_types": ["Tissue-resident T cell", "NKT cell"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "C-X-C motif chemokine receptor 6. Marks tissue-resident memory T cells and NKT cells. Mediates liver homing and retention in sinusoids."},
    # --- Proliferation markers ---
    {"gene": "MKI67", "cell_types": ["Proliferating cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Marker of proliferation Ki-67. Nuclear protein present in all active cell cycle phases (G1/S/G2/M) but absent in G0. Gold standard proliferation marker."},
    {"gene": "TOP2A", "cell_types": ["Proliferating cell"], "specificity": "high", "source_db": "panglaodb",
     "description": "DNA topoisomerase II alpha. Essential enzyme for DNA replication and chromosome segregation. Marks S/G2/M phase cells."},
    {"gene": "CDK1", "cell_types": ["Proliferating cell"], "specificity": "high", "source_db": "panglaodb",
     "description": "Cyclin-dependent kinase 1. Master regulator of M-phase entry. Required for mitosis in all cell types."},
    {"gene": "PCNA", "cell_types": ["Proliferating cell"], "specificity": "moderate", "source_db": "cellmarker",
     "description": "Proliferating cell nuclear antigen. DNA sliding clamp for replicative DNA polymerases. Expressed in S phase and DNA repair."},
    {"gene": "MCM2", "cell_types": ["Proliferating cell"], "specificity": "moderate", "source_db": "panglaodb",
     "description": "Minichromosome maintenance complex component 2. Essential DNA replication licensing factor marking cells that have exited G0."},
    # --- Mast cell markers ---
    {"gene": "TPSB2", "cell_types": ["Mast cell"], "specificity": "high", "source_db": "cellmarker",
     "description": "Tryptase beta 2. Serine protease stored in mast cell secretory granules. Distinguishes connective tissue mast cells (TPSB2+TPSAB1+) from mucosal mast cells."},
]


def get_marker_count() -> int:
    """Return the number of seed marker gene records."""
    return len(MARKER_GENE_RECORDS)


def get_marker_sources() -> List[str]:
    """Return the unique source databases from seed data."""
    sources = sorted(set(r["source_db"] for r in MARKER_GENE_RECORDS))
    return sources


# ===================================================================
# MARKER PARSER
# ===================================================================


class MarkerParser(BaseIngestParser):
    """Ingest parser for CellMarker / PanglaoDB marker gene data.

    Seeds the knowledge base with canonical marker genes, their cell type
    associations, specificity ratings, and functional descriptions.

    Usage::

        parser = MarkerParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        collection_manager=None,
        embedder=None,
    ) -> None:
        super().__init__(
            source_name="cellmarker",
            collection_manager=collection_manager,
            embedder=embedder,
        )

    def fetch(self, **kwargs) -> List[Dict[str, Any]]:
        """Return seed marker gene records.

        In production this would query CellMarker2.0 and PanglaoDB APIs.
        For seed mode it returns curated records.
        """
        return list(MARKER_GENE_RECORDS)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse marker gene dictionaries into IngestRecord objects."""
        records = []
        for entry in raw_data:
            cell_types_str = ", ".join(entry.get("cell_types", []))
            text = (
                f"Marker gene: {entry['gene']}. "
                f"Cell types: {cell_types_str}. "
                f"Specificity: {entry.get('specificity', 'unknown')}. "
                f"Source: {entry.get('source_db', 'unknown')}. "
                f"{entry.get('description', '')}"
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "gene": entry.get("gene", ""),
                    "cell_types": entry.get("cell_types", []),
                    "specificity": entry.get("specificity", "unknown"),
                    "source_db": entry.get("source_db", ""),
                },
                collection_name="sc_markers",
                record_id=f"marker_{entry.get('gene', 'unknown')}",
                source=entry.get("source_db", "cellmarker"),
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate that a marker record has minimum required data."""
        if len(record.text) < 20:
            return False
        if not record.metadata.get("gene"):
            return False
        if not record.metadata.get("cell_types"):
            return False
        return True
