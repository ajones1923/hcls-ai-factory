"""CellxGene / Human Cell Atlas cell type parser for the Single-Cell Intelligence Agent.

Seeds 49 cell type records from the CellxGene Census and Human Cell Atlas
reference datasets, covering major cell lineages, tissue-resident populations,
immune subtypes, and stem/progenitor cells.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: 49 CELL TYPE RECORDS FROM CellxGene / HCA
# ===================================================================

CELL_TYPE_RECORDS: List[Dict[str, Any]] = [
    # --- Immune: T cells ---
    {
        "cell_type_id": "CL:0000084",
        "cell_type": "T cell",
        "lineage": "lymphoid",
        "tissue": "blood, lymph node, spleen, thymus",
        "markers": ["CD3D", "CD3E", "CD3G", "TRAC"],
        "description": "T lymphocytes are adaptive immune cells that mediate cellular immunity. "
                       "They recognize antigens via T cell receptor (TCR) complexes and differentiate "
                       "into CD4+ helper, CD8+ cytotoxic, and regulatory T cell subsets.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000624",
        "cell_type": "CD4-positive helper T cell",
        "lineage": "lymphoid",
        "tissue": "blood, lymph node, spleen",
        "markers": ["CD4", "IL7R", "TCF7", "LEF1", "CCR7"],
        "description": "CD4+ T helper cells coordinate adaptive immune responses by producing "
                       "cytokines that activate B cells, macrophages, and cytotoxic T cells. "
                       "Subtypes include Th1, Th2, Th17, and Tfh cells.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000625",
        "cell_type": "CD8-positive cytotoxic T cell",
        "lineage": "lymphoid",
        "tissue": "blood, lymph node, tumor",
        "markers": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7"],
        "description": "CD8+ cytotoxic T lymphocytes directly kill virus-infected and tumor cells "
                       "via perforin-granzyme and Fas/FasL pathways. Key effectors in anti-tumor immunity.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000815",
        "cell_type": "Regulatory T cell",
        "lineage": "lymphoid",
        "tissue": "blood, lymph node, tumor",
        "markers": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
        "description": "Regulatory T cells (Tregs) suppress immune responses and maintain "
                       "self-tolerance. FOXP3 is the master transcription factor. Tregs are enriched "
                       "in the tumor microenvironment and can suppress anti-tumor immunity.",
        "source": "cellxgene",
    },
    # --- Immune: B cells ---
    {
        "cell_type_id": "CL:0000236",
        "cell_type": "B cell",
        "lineage": "lymphoid",
        "tissue": "blood, bone marrow, lymph node, spleen",
        "markers": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"],
        "description": "B lymphocytes mediate humoral immunity by producing antibodies. "
                       "They undergo class switching and somatic hypermutation in germinal centers.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000786",
        "cell_type": "Plasma cell",
        "lineage": "lymphoid",
        "tissue": "bone marrow, lymph node, mucosa",
        "markers": ["SDC1", "TNFRSF17", "XBP1", "PRDM1", "IRF4"],
        "description": "Terminally differentiated B cells that are the primary antibody-secreting "
                       "cells. Long-lived plasma cells reside in bone marrow niches.",
        "source": "cellxgene",
    },
    # --- Immune: NK cells ---
    {
        "cell_type_id": "CL:0000623",
        "cell_type": "Natural killer cell",
        "lineage": "lymphoid",
        "tissue": "blood, liver, uterus, spleen",
        "markers": ["NCAM1", "NKG7", "KLRD1", "GNLY", "FCGR3A"],
        "description": "NK cells are innate lymphoid cells that kill virus-infected and tumor cells "
                       "without prior antigen sensitization. CD56bright and CD56dim subsets differ "
                       "in cytokine production and cytotoxicity.",
        "source": "cellxgene",
    },
    # --- Myeloid ---
    {
        "cell_type_id": "CL:0000235",
        "cell_type": "Macrophage",
        "lineage": "myeloid",
        "tissue": "lung, liver, brain, adipose, peritoneum",
        "markers": ["CD68", "CD163", "MRC1", "CSF1R", "MARCO"],
        "description": "Macrophages are tissue-resident phagocytes that clear debris, present "
                       "antigens, and regulate inflammation. M1 (pro-inflammatory) and M2 "
                       "(anti-inflammatory) polarization states influence tissue repair and tumor progression.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000451",
        "cell_type": "Dendritic cell",
        "lineage": "myeloid",
        "tissue": "blood, skin, lymph node, lung",
        "markers": ["ITGAX", "HLA-DRA", "FLT3", "CLEC9A", "CD1C"],
        "description": "Professional antigen-presenting cells that bridge innate and adaptive "
                       "immunity. Conventional DC1 (cDC1) cross-present antigens to CD8+ T cells; "
                       "cDC2 activate CD4+ T cells.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000576",
        "cell_type": "Monocyte",
        "lineage": "myeloid",
        "tissue": "blood, bone marrow",
        "markers": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN"],
        "description": "Circulating monocytes differentiate into macrophages and dendritic cells "
                       "upon tissue infiltration. Classical (CD14++CD16-), intermediate, and "
                       "non-classical subtypes have distinct functions.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000775",
        "cell_type": "Neutrophil",
        "lineage": "myeloid",
        "tissue": "blood, bone marrow, inflamed tissue",
        "markers": ["CSF3R", "FCGR3B", "CXCR2", "S100A8", "MMP9"],
        "description": "Most abundant circulating leukocyte. First responders to infection via "
                       "phagocytosis, degranulation, and neutrophil extracellular trap (NET) formation. "
                       "Tumor-associated neutrophils can be N1 (anti-tumor) or N2 (pro-tumor).",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000767",
        "cell_type": "Mast cell",
        "lineage": "myeloid",
        "tissue": "skin, gut, lung, connective tissue",
        "markers": ["KIT", "TPSAB1", "CPA3", "FCER1A", "HDC"],
        "description": "Tissue-resident granulocytes containing histamine and proteases. "
                       "Central to allergic responses and IgE-mediated immunity. Also involved "
                       "in tumor angiogenesis.",
        "source": "cellxgene",
    },
    # --- Epithelial ---
    {
        "cell_type_id": "CL:0000066",
        "cell_type": "Epithelial cell",
        "lineage": "epithelial",
        "tissue": "lung, intestine, skin, kidney, breast",
        "markers": ["EPCAM", "KRT18", "KRT8", "CDH1", "CLDN4"],
        "description": "Epithelial cells line body surfaces and cavities. They provide barrier "
                       "function, secrete mucus and surfactant, and participate in absorption. "
                       "Many carcinomas arise from epithelial cells.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0002633",
        "cell_type": "Alveolar type II cell",
        "lineage": "epithelial",
        "tissue": "lung",
        "markers": ["SFTPC", "SFTPB", "ABCA3", "NKX2-1", "LAMP3"],
        "description": "AT2 cells produce pulmonary surfactant and serve as progenitors for "
                       "alveolar type I cells during lung repair. Key cells in COVID-19 pathology "
                       "and lung adenocarcinoma.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0002563",
        "cell_type": "Intestinal epithelial cell",
        "lineage": "epithelial",
        "tissue": "small intestine, colon",
        "markers": ["EPCAM", "VIL1", "FABP2", "CDX2", "MUC2"],
        "description": "Intestinal epithelial cells include absorptive enterocytes, goblet cells, "
                       "Paneth cells, enteroendocrine cells, and Lgr5+ stem cells organized in "
                       "crypt-villus architecture.",
        "source": "cellxgene",
    },
    # --- Stromal ---
    {
        "cell_type_id": "CL:0000057",
        "cell_type": "Fibroblast",
        "lineage": "mesenchymal",
        "tissue": "skin, lung, heart, tumor stroma",
        "markers": ["COL1A1", "DCN", "LUM", "PDGFRA", "FAP"],
        "description": "Fibroblasts produce extracellular matrix and growth factors. "
                       "Cancer-associated fibroblasts (CAFs) remodel the tumor stroma and can "
                       "promote or restrain tumor growth depending on subtype.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0000669",
        "cell_type": "Pericyte",
        "lineage": "mesenchymal",
        "tissue": "vasculature, brain, retina",
        "markers": ["PDGFRB", "RGS5", "CSPG4", "ACTA2", "NOTCH3"],
        "description": "Mural cells that wrap around endothelial cells in capillaries and "
                       "venules. They regulate blood flow, vessel stability, and blood-brain "
                       "barrier integrity.",
        "source": "cellxgene",
    },
    # --- Endothelial ---
    {
        "cell_type_id": "CL:0000115",
        "cell_type": "Endothelial cell",
        "lineage": "endothelial",
        "tissue": "vasculature, liver, lung, brain",
        "markers": ["PECAM1", "VWF", "CDH5", "ERG", "KDR"],
        "description": "Endothelial cells line blood and lymphatic vessels. They regulate "
                       "vascular tone, permeability, and leukocyte trafficking. Tumor endothelium "
                       "is phenotypically distinct from normal vasculature.",
        "source": "cellxgene",
    },
    # --- Stem / Progenitor ---
    {
        "cell_type_id": "CL:0000037",
        "cell_type": "Hematopoietic stem cell",
        "lineage": "stem",
        "tissue": "bone marrow, cord blood",
        "markers": ["CD34", "KIT", "THY1", "CXCR4", "PTPRC"],
        "description": "Self-renewing multipotent cells that give rise to all blood cell lineages. "
                       "HSCs reside in bone marrow niches and are used in transplantation for "
                       "hematologic malignancies.",
        "source": "cellxgene",
    },
    {
        "cell_type_id": "CL:0002322",
        "cell_type": "Embryonic stem cell",
        "lineage": "stem",
        "tissue": "inner cell mass",
        "markers": ["POU5F1", "NANOG", "SOX2", "LIN28A", "DPPA4"],
        "description": "Pluripotent stem cells derived from the inner cell mass of blastocysts. "
                       "They can differentiate into all three germ layers and are widely used in "
                       "regenerative medicine research.",
        "source": "hca",
    },
    # --- Neural ---
    {
        "cell_type_id": "CL:0000540",
        "cell_type": "Neuron",
        "lineage": "neural",
        "tissue": "brain, spinal cord, peripheral nerve",
        "markers": ["RBFOX3", "SYP", "MAP2", "SNAP25", "SLC17A7"],
        "description": "Electrically excitable cells that transmit signals via synapses. "
                       "Subtypes include glutamatergic, GABAergic, dopaminergic, and "
                       "serotonergic neurons.",
        "source": "hca",
    },
    {
        "cell_type_id": "CL:0000127",
        "cell_type": "Astrocyte",
        "lineage": "neural",
        "tissue": "brain, spinal cord",
        "markers": ["GFAP", "AQP4", "SLC1A3", "S100B", "ALDH1L1"],
        "description": "Star-shaped glial cells that support neurons, maintain the blood-brain "
                       "barrier, regulate neurotransmitter levels, and participate in "
                       "neuroinflammation.",
        "source": "hca",
    },
    {
        "cell_type_id": "CL:0000129",
        "cell_type": "Microglial cell",
        "lineage": "myeloid",
        "tissue": "brain",
        "markers": ["CX3CR1", "TMEM119", "P2RY12", "ITGAM", "AIF1"],
        "description": "Resident macrophages of the central nervous system derived from yolk sac "
                       "progenitors. They survey the brain parenchyma, prune synapses during "
                       "development, and respond to neurodegeneration.",
        "source": "cellxgene",
    },
    # --- Muscle ---
    {
        "cell_type_id": "CL:0000187",
        "cell_type": "Muscle cell",
        "lineage": "mesenchymal",
        "tissue": "skeletal muscle, heart, smooth muscle",
        "markers": ["DES", "ACTA1", "MYH11", "TTN", "TNNT2"],
        "description": "Contractile cells including skeletal myocytes, cardiomyocytes, and "
                       "smooth muscle cells. Skeletal muscle satellite cells are tissue-resident "
                       "stem cells for muscle regeneration.",
        "source": "cellxgene",
    },
    # --- Hepatocyte ---
    {
        "cell_type_id": "CL:0000182",
        "cell_type": "Hepatocyte",
        "lineage": "epithelial",
        "tissue": "liver",
        "markers": ["ALB", "HNF4A", "SERPINA1", "CYP3A4", "APOB"],
        "description": "Primary parenchymal cells of the liver responsible for metabolism, "
                       "detoxification, protein synthesis, and bile production. Periportal and "
                       "pericentral zones have distinct metabolic programs.",
        "source": "cellxgene",
    },
    # --- Adipocyte ---
    {
        "cell_type_id": "CL:0000136",
        "cell_type": "Adipocyte",
        "lineage": "mesenchymal",
        "tissue": "adipose tissue",
        "markers": ["ADIPOQ", "LEP", "FABP4", "PPARG", "PLIN1"],
        "description": "Lipid-storing cells that regulate energy homeostasis and secrete "
                       "adipokines. White and brown adipocyte subtypes have distinct metabolic roles.",
        "source": "cellxgene",
    },
    # --- Platelet ---
    {
        "cell_type_id": "CL:0000233",
        "cell_type": "Platelet",
        "lineage": "myeloid",
        "tissue": "blood",
        "markers": ["ITGA2B", "GP1BA", "PF4", "PPBP", "SELP"],
        "description": "Anucleate cell fragments derived from megakaryocytes. They mediate "
                       "hemostasis, wound healing, and can interact with tumor cells to promote "
                       "metastasis.",
        "source": "cellxgene",
    },
    # --- Erythrocyte precursor ---
    {
        "cell_type_id": "CL:0000764",
        "cell_type": "Erythroid progenitor",
        "lineage": "erythroid",
        "tissue": "bone marrow, fetal liver",
        "markers": ["GYPA", "KLF1", "HBB", "TFRC", "GATA1"],
        "description": "Erythroid progenitor cells differentiate through proerythroblast, "
                       "basophilic, polychromatic, and orthochromatic stages to produce reticulocytes "
                       "and mature red blood cells.",
        "source": "cellxgene",
    },
    # --- Melanocyte ---
    {
        "cell_type_id": "CL:0000148",
        "cell_type": "Melanocyte",
        "lineage": "neural_crest",
        "tissue": "skin, eye, hair follicle",
        "markers": ["MITF", "TYR", "PMEL", "DCT", "MLANA"],
        "description": "Neural crest-derived pigment-producing cells. Melanoma arises from "
                       "malignant transformation of melanocytes, with BRAF V600E and NRAS mutations "
                       "as common drivers.",
        "source": "cellxgene",
    },
    # --- Chondrocyte ---
    {
        "cell_type_id": "CL:0000138",
        "cell_type": "Chondrocyte",
        "lineage": "mesenchymal",
        "tissue": "cartilage, growth plate",
        "markers": ["SOX9", "COL2A1", "ACAN", "COL11A1", "COMP"],
        "description": "Cartilage-resident cells that produce and maintain the cartilaginous "
                       "matrix. Important in joint disease, skeletal development, and tissue "
                       "engineering applications.",
        "source": "cellxgene",
    },
    # --- Osteoblast ---
    {
        "cell_type_id": "CL:0000062",
        "cell_type": "Osteoblast",
        "lineage": "mesenchymal",
        "tissue": "bone",
        "markers": ["RUNX2", "SP7", "BGLAP", "COL1A1", "ALPL"],
        "description": "Bone-forming cells derived from mesenchymal stem cells. They synthesize "
                       "osteoid matrix and regulate bone mineralization.",
        "source": "cellxgene",
    },
    # --- Mesenchymal stem cell ---
    {
        "cell_type_id": "CL:0000134",
        "cell_type": "Mesenchymal stem cell",
        "lineage": "mesenchymal",
        "tissue": "bone marrow, adipose, umbilical cord",
        "markers": ["ENG", "NT5E", "THY1", "VCAM1", "LEPR"],
        "description": "Multipotent stromal cells that differentiate into osteoblasts, "
                       "chondrocytes, and adipocytes. Used in cell therapy and tissue engineering.",
        "source": "cellxgene",
    },
    # --- Trophoblast ---
    {
        "cell_type_id": "CL:0000351",
        "cell_type": "Trophoblast",
        "lineage": "extraembryonic",
        "tissue": "placenta",
        "markers": ["KRT7", "GATA3", "TFAP2A", "HLA-G", "CGA"],
        "description": "Cells forming the outer layer of the blastocyst that mediate "
                       "implantation and placental development. Subtypes include cytotrophoblasts, "
                       "syncytiotrophoblasts, and extravillous trophoblasts.",
        "source": "hca",
    },
    # --- Oligodendrocyte ---
    {
        "cell_type_id": "CL:0000128",
        "cell_type": "Oligodendrocyte",
        "lineage": "neural",
        "tissue": "brain, spinal cord",
        "markers": ["MBP", "PLP1", "MOG", "MAG", "OLIG2"],
        "description": "Myelinating glial cells of the central nervous system. Each oligodendrocyte "
                       "can myelinate multiple axon segments. Loss of oligodendrocytes contributes to "
                       "demyelinating diseases like multiple sclerosis.",
        "source": "hca",
    },
    # --- Plasmacytoid dendritic cell ---
    {
        "cell_type_id": "CL:0000784",
        "cell_type": "Plasmacytoid dendritic cell",
        "lineage": "myeloid",
        "tissue": "blood, lymphoid tissue",
        "markers": ["CLEC4C", "IL3RA", "IRF7", "TCF4", "LILRA4"],
        "description": "Plasmacytoid dendritic cells (pDCs) are specialized innate immune cells that "
                       "produce large amounts of type I interferon (IFN-alpha/beta) in response to viral "
                       "nucleic acids. They express TLR7 and TLR9 and bridge innate antiviral defense "
                       "with adaptive immunity.",
        "source": "cellxgene",
    },
    # --- cDC1 ---
    {
        "cell_type_id": "CL:0002399",
        "cell_type": "Conventional dendritic cell type 1 (cDC1)",
        "lineage": "myeloid",
        "tissue": "blood, lymph node, tumor",
        "markers": ["CLEC9A", "XCR1", "BATF3", "IRF8", "THBD"],
        "description": "cDC1s are specialized in cross-presentation of exogenous antigens on MHC class I "
                       "to CD8+ T cells. They are critical for anti-tumor immunity and are regulated by "
                       "the transcription factors BATF3 and IRF8. Their presence in tumors correlates "
                       "with immunotherapy response.",
        "source": "cellxgene",
    },
    # --- cDC2 ---
    {
        "cell_type_id": "CL:0002394",
        "cell_type": "Conventional dendritic cell type 2 (cDC2)",
        "lineage": "myeloid",
        "tissue": "blood, lymph node, skin, lung",
        "markers": ["CD1C", "FCER1A", "CLEC10A", "IRF4", "SIRPA"],
        "description": "cDC2s present antigens on MHC class II to activate CD4+ T helper cells. They "
                       "are more heterogeneous than cDC1s and can polarize Th1, Th2, and Th17 responses "
                       "depending on the inflammatory context.",
        "source": "cellxgene",
    },
    # --- ILC1 ---
    {
        "cell_type_id": "CL:0001069",
        "cell_type": "Innate lymphoid cell type 1 (ILC1)",
        "lineage": "lymphoid",
        "tissue": "gut, liver, tonsil",
        "markers": ["TBX21", "IFNG", "IL12RB2", "NCR1", "EOMES"],
        "description": "ILC1s are tissue-resident innate lymphoid cells that produce IFN-gamma and TNF "
                       "in response to intracellular pathogens. They are T-bet dependent and contribute "
                       "to mucosal defense in the gut and liver.",
        "source": "cellxgene",
    },
    # --- ILC2 ---
    {
        "cell_type_id": "CL:0001070",
        "cell_type": "Innate lymphoid cell type 2 (ILC2)",
        "lineage": "lymphoid",
        "tissue": "lung, skin, gut, adipose",
        "markers": ["GATA3", "IL5", "IL13", "KLRG1", "AREG"],
        "description": "ILC2s are tissue-resident innate lymphoid cells that produce type 2 cytokines "
                       "(IL-5, IL-13) in response to epithelial alarmins (IL-25, IL-33, TSLP). They "
                       "drive allergic inflammation and participate in tissue repair and metabolic "
                       "homeostasis.",
        "source": "cellxgene",
    },
    # --- ILC3 ---
    {
        "cell_type_id": "CL:0001071",
        "cell_type": "Innate lymphoid cell type 3 (ILC3)",
        "lineage": "lymphoid",
        "tissue": "gut, tonsil, lymph node",
        "markers": ["RORC", "IL17A", "IL22", "NCR2", "KIT"],
        "description": "ILC3s are RORgamma-t-dependent innate lymphoid cells that produce IL-17A and "
                       "IL-22 to maintain mucosal barrier integrity. They regulate intestinal homeostasis, "
                       "lymphoid tissue organogenesis, and defense against extracellular bacteria and fungi.",
        "source": "cellxgene",
    },
    # --- MAIT cell ---
    {
        "cell_type_id": "CL:0000940",
        "cell_type": "Mucosal-associated invariant T cell (MAIT)",
        "lineage": "lymphoid",
        "tissue": "liver, blood, gut, lung",
        "markers": ["TRAV1-2", "SLC4A10", "KLRB1", "DPP4", "IL18R1"],
        "description": "MAIT cells are innate-like T cells that express a semi-invariant TCR alpha chain "
                       "(TRAV1-2/TRAJ33) and recognize microbial riboflavin metabolites presented by MR1. "
                       "They are enriched in the liver and mucosal tissues and respond rapidly to bacterial "
                       "infection.",
        "source": "cellxgene",
    },
    # --- Gamma-delta T cell ---
    {
        "cell_type_id": "CL:0000798",
        "cell_type": "Gamma-delta T cell",
        "lineage": "lymphoid",
        "tissue": "skin, gut, blood, reproductive tract",
        "markers": ["TRGV9", "TRDV2", "TRDC", "CD3E", "IL17A"],
        "description": "Gamma-delta T cells express a TCR composed of gamma and delta chains and bridge "
                       "innate and adaptive immunity. Vgamma9Vdelta2 T cells recognize phosphoantigens "
                       "from microbes and stressed cells. They are enriched in epithelial tissues and "
                       "contribute to anti-tumor surveillance.",
        "source": "cellxgene",
    },
    # --- Mast cell (additional ontology entry) ---
    {
        "cell_type_id": "CL:0000097",
        "cell_type": "Mast cell (connective tissue type)",
        "lineage": "myeloid",
        "tissue": "connective tissue, skin, peritoneum",
        "markers": ["KIT", "FCER1A", "TPSB2", "CPA3", "HPGDS"],
        "description": "Connective tissue mast cells (CTMC) are tissue-resident granulocytes that "
                       "contain tryptase and chymase. They are found in skin, peritoneum, and connective "
                       "tissue and play roles in allergic reactions, fibrosis, and wound healing.",
        "source": "cellxgene",
    },
    # --- Megakaryocyte ---
    {
        "cell_type_id": "CL:0000556",
        "cell_type": "Megakaryocyte",
        "lineage": "myeloid",
        "tissue": "bone marrow, lung",
        "markers": ["ITGA2B", "GP1BA", "PF4", "PPBP", "GATA1"],
        "description": "Megakaryocytes are large polyploid cells in the bone marrow that produce "
                       "platelets through cytoplasmic fragmentation. They undergo endomitosis to "
                       "achieve polyploidy and extend proplatelets into sinusoidal blood vessels. "
                       "Thrombopoietin (TPO) is the primary growth factor for megakaryopoiesis.",
        "source": "cellxgene",
    },
    # --- Erythroid progenitor (early) ---
    {
        "cell_type_id": "CL:0000038",
        "cell_type": "Erythroid progenitor (early)",
        "lineage": "erythroid",
        "tissue": "bone marrow, fetal liver",
        "markers": ["GYPA", "KLF1", "GATA1", "TAL1", "EPOR"],
        "description": "Early erythroid progenitors (BFU-E and CFU-E) are committed progenitors that "
                       "give rise to erythroblasts under the influence of erythropoietin. They undergo "
                       "chromatin condensation and enucleation during terminal erythroid differentiation.",
        "source": "cellxgene",
    },
    # --- Goblet cell ---
    {
        "cell_type_id": "CL:0000160",
        "cell_type": "Goblet cell",
        "lineage": "epithelial",
        "tissue": "intestine, airway, conjunctiva",
        "markers": ["MUC2", "TFF3", "SPDEF", "FCGBP", "CLCA1"],
        "description": "Goblet cells are specialized secretory epithelial cells that produce and secrete "
                       "mucins to form the protective mucus barrier in the intestine and airways. They "
                       "differentiate from Lgr5+ stem cells via Notch pathway inhibition and are "
                       "regulated by the transcription factor SPDEF.",
        "source": "cellxgene",
    },
    # --- Schwann cell ---
    {
        "cell_type_id": "CL:0002573",
        "cell_type": "Schwann cell",
        "lineage": "neural_crest",
        "tissue": "peripheral nerve",
        "markers": ["MPZ", "PMP22", "SOX10", "S100B", "EGR2"],
        "description": "Schwann cells are the principal glial cells of the peripheral nervous system. "
                       "Myelinating Schwann cells wrap individual axons with myelin sheaths, while "
                       "non-myelinating Schwann cells (Remak cells) ensheathe multiple small-caliber axons. "
                       "They are essential for nerve regeneration after injury.",
        "source": "cellxgene",
    },
    # --- Podocyte ---
    {
        "cell_type_id": "CL:0000653",
        "cell_type": "Podocyte",
        "lineage": "epithelial",
        "tissue": "kidney glomerulus",
        "markers": ["NPHS1", "NPHS2", "WT1", "PODXL", "SYNPO"],
        "description": "Podocytes are highly specialized epithelial cells in the kidney glomerulus that "
                       "form the final barrier of the glomerular filtration apparatus. They extend "
                       "interdigitating foot processes connected by slit diaphragms. Podocyte injury "
                       "or loss is central to proteinuria and glomerular disease.",
        "source": "cellxgene",
    },
]


def get_cell_type_count() -> int:
    """Return the number of seed cell type records."""
    return len(CELL_TYPE_RECORDS)


def get_cell_lineages() -> List[str]:
    """Return the unique lineage categories from seed data."""
    lineages = sorted(set(r["lineage"] for r in CELL_TYPE_RECORDS))
    return lineages


# ===================================================================
# CellxGene PARSER
# ===================================================================


class CellxGeneParser(BaseIngestParser):
    """Ingest parser for CellxGene / Human Cell Atlas cell type data.

    Seeds the knowledge base with curated cell type definitions,
    canonical markers, tissue distributions, and functional descriptions.

    Usage::

        parser = CellxGeneParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        collection_manager=None,
        embedder=None,
    ) -> None:
        super().__init__(
            source_name="cellxgene",
            collection_manager=collection_manager,
            embedder=embedder,
        )

    def fetch(self, **kwargs) -> List[Dict[str, Any]]:
        """Return seed cell type records.

        In production this would query the CellxGene Census API.
        For seed mode it returns curated records.
        """
        return list(CELL_TYPE_RECORDS)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse cell type dictionaries into IngestRecord objects."""
        records = []
        for entry in raw_data:
            markers_str = ", ".join(entry.get("markers", []))
            text = (
                f"Cell type: {entry['cell_type']} ({entry.get('cell_type_id', '')}). "
                f"Lineage: {entry.get('lineage', '')}. "
                f"Tissue: {entry.get('tissue', '')}. "
                f"Canonical markers: {markers_str}. "
                f"{entry.get('description', '')}"
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "cell_type_id": entry.get("cell_type_id", ""),
                    "cell_type": entry.get("cell_type", ""),
                    "lineage": entry.get("lineage", ""),
                    "tissue": entry.get("tissue", ""),
                    "markers": entry.get("markers", []),
                },
                collection_name="sc_cell_types",
                record_id=entry.get("cell_type_id", ""),
                source=entry.get("source", "cellxgene"),
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate that a cell type record has minimum required data."""
        if len(record.text) < 30:
            return False
        if not record.metadata.get("cell_type"):
            return False
        if not record.metadata.get("markers"):
            return False
        return True
