"""Single-Cell Intelligence Agent -- Query Expansion Module.

Provides entity alias resolution, synonym mapping, and query expansion
for single-cell genomics queries.  Ensures that abbreviations, technology
names, cell type aliases, and colloquial terms are normalized and expanded
to improve vector-search recall in the RAG pipeline.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


# ===================================================================
# ENTITY ALIASES (150+ single-cell-specific aliases)
# ===================================================================

ENTITY_ALIASES: Dict[str, str] = {
    # -- Technology abbreviations --
    "scRNA-seq": "single-cell RNA sequencing",
    "scRNAseq": "single-cell RNA sequencing",
    "snRNA-seq": "single-nucleus RNA sequencing",
    "snRNAseq": "single-nucleus RNA sequencing",
    "scATAC-seq": "single-cell ATAC sequencing",
    "scATACseq": "single-cell ATAC sequencing",
    "CITE-seq": "Cellular Indexing of Transcriptomes and Epitopes by Sequencing",
    "CITEseq": "Cellular Indexing of Transcriptomes and Epitopes by Sequencing",
    "ECCITE-seq": "expanded CITE-seq",
    "SHARE-seq": "simultaneous high-throughput ATAC and RNA expression with sequencing",
    "TEA-seq": "transcription, epitopes, and accessibility sequencing",
    "10x": "10x Genomics Chromium",
    "10X": "10x Genomics Chromium",
    "GEX": "gene expression",
    "ADT": "antibody-derived tag",
    "HTO": "hashtag oligonucleotide",
    "VDJ": "variable diversity joining",
    "TCR-seq": "T-cell receptor sequencing",
    "BCR-seq": "B-cell receptor sequencing",
    "CROP-seq": "CRISPR droplet sequencing",
    "Perturb-seq": "perturbation with single-cell sequencing",
    "sci-RNA-seq": "single-cell combinatorial indexing RNA sequencing",
    "Smart-seq2": "switching mechanism at 5-prime end of RNA template version 2",
    "CEL-Seq2": "cell expression by linear amplification and sequencing 2",
    "Drop-seq": "droplet-based single-cell sequencing",
    "inDrop": "indexing droplets single-cell sequencing",
    "MARS-seq": "massively parallel single-cell RNA sequencing",
    "Seq-Well": "portable single-cell sequencing platform",
    "Slide-seq": "slide-based spatial transcriptomics",
    "STARmap": "spatially resolved transcript amplification readout mapping",
    "osmFISH": "ouroboros single molecule FISH",
    "seqFISH": "sequential fluorescence in situ hybridization",
    "seqFISH+": "enhanced sequential FISH",

    # -- Spatial platforms --
    "MERFISH": "multiplexed error-robust fluorescence in situ hybridization",
    "MERSCOPE": "Vizgen MERSCOPE spatial genomics platform",
    "Visium": "10x Genomics Visium spatial transcriptomics",
    "Visium HD": "10x Genomics Visium HD single-cell resolution spatial",
    "CytAssist": "10x Genomics CytAssist FFPE spatial",
    "Xenium": "10x Genomics Xenium in situ platform",
    "CODEX": "co-detection by indexing spatial proteomics",
    "CosMx": "NanoString CosMx spatial molecular imager",
    "Stereo-seq": "BGI spatiotemporal enhanced resolution omics sequencing",
    "DBiT-seq": "deterministic barcoding in tissue for spatial omics",
    "IMC": "imaging mass cytometry",
    "MIBI": "multiplexed ion beam imaging",
    "MIBI-TOF": "multiplexed ion beam imaging by time of flight",

    # -- Cell type abbreviations --
    "Treg": "regulatory T cell",
    "Tregs": "regulatory T cells",
    "CTL": "cytotoxic T lymphocyte",
    "CTLs": "cytotoxic T lymphocytes",
    "TIL": "tumor-infiltrating lymphocyte",
    "TILs": "tumor-infiltrating lymphocytes",
    "NKT": "natural killer T cell",
    "MAIT": "mucosal-associated invariant T cell",
    "Tfh": "T follicular helper cell",
    "Th1": "T helper 1 cell",
    "Th2": "T helper 2 cell",
    "Th17": "T helper 17 cell",
    "cDC1": "conventional dendritic cell type 1",
    "cDC2": "conventional dendritic cell type 2",
    "pDC": "plasmacytoid dendritic cell",
    "moDC": "monocyte-derived dendritic cell",
    "TAM": "tumor-associated macrophage",
    "TAMs": "tumor-associated macrophages",
    "CAF": "cancer-associated fibroblast",
    "CAFs": "cancer-associated fibroblasts",
    "MDSC": "myeloid-derived suppressor cell",
    "MDSCs": "myeloid-derived suppressor cells",
    "HSC": "hematopoietic stem cell",
    "HSCs": "hematopoietic stem cells",
    "MSC": "mesenchymal stem cell",
    "MSCs": "mesenchymal stem cells",
    "OPC": "oligodendrocyte precursor cell",
    "OPCs": "oligodendrocyte precursor cells",
    "AT1": "alveolar type 1 cell",
    "AT2": "alveolar type 2 cell",
    "GC B": "germinal center B cell",
    "ILC": "innate lymphoid cell",
    "ILCs": "innate lymphoid cells",
    "ILC1": "innate lymphoid cell type 1",
    "ILC2": "innate lymphoid cell type 2",
    "ILC3": "innate lymphoid cell type 3",
    "Trm": "tissue-resident memory T cell",
    "PBMC": "peripheral blood mononuclear cell",
    "PBMCs": "peripheral blood mononuclear cells",

    # -- Gene / marker abbreviations --
    "PD-1": "PDCD1 programmed death 1",
    "PD1": "PDCD1 programmed death 1",
    "PD-L1": "CD274 programmed death ligand 1",
    "PDL1": "CD274 programmed death ligand 1",
    "CTLA-4": "CTLA4 cytotoxic T-lymphocyte associated protein 4",
    "TIM-3": "HAVCR2 T-cell immunoglobulin mucin 3",
    "LAG-3": "LAG3 lymphocyte activation gene 3",
    "SMA": "ACTA2 alpha-smooth muscle actin",
    "aSMA": "ACTA2 alpha-smooth muscle actin",
    "NeuN": "RBFOX3 neuronal nuclei antigen",
    "CD56": "NCAM1 neural cell adhesion molecule",
    "CD138": "SDC1 syndecan 1",
    "CD31": "PECAM1 platelet endothelial cell adhesion molecule",
    "CD20": "MS4A1 membrane spanning 4-domains A1",
    "CD25": "IL2RA interleukin 2 receptor subunit alpha",
    "CD117": "KIT stem cell factor receptor",
    "CD45": "PTPRC protein tyrosine phosphatase receptor type C",
    "EpCAM": "EPCAM epithelial cell adhesion molecule",

    # -- Analysis method abbreviations --
    "UMAP": "uniform manifold approximation and projection",
    "tSNE": "t-distributed stochastic neighbor embedding",
    "t-SNE": "t-distributed stochastic neighbor embedding",
    "PCA": "principal component analysis",
    "HVG": "highly variable genes",
    "HVGs": "highly variable genes",
    "DEG": "differentially expressed gene",
    "DEGs": "differentially expressed genes",
    "DE": "differential expression",
    "DGE": "differential gene expression",
    "GRN": "gene regulatory network",
    "GRNs": "gene regulatory networks",
    "CNV": "copy number variation",
    "CNVs": "copy number variations",
    "CNA": "copy number alteration",
    "QC": "quality control",
    "UMI": "unique molecular identifier",
    "UMIs": "unique molecular identifiers",
    "SCT": "sctransform normalization",
    "RNA velocity": "RNA velocity spliced/unspliced ratio analysis",
    "PAGA": "partition-based graph abstraction",
    "SCENIC": "single-cell regulatory network inference and clustering",
    "SCENIC+": "enhanced SCENIC with chromatin accessibility",
    "CellChat": "cell-cell communication analysis tool",
    "CellPhoneDB": "cell-cell interaction database",
    "NicheNet": "ligand-receptor activity modeling",
    "Monocle": "pseudotime trajectory analysis",
    "Monocle3": "Monocle version 3 trajectory analysis",
    "Palantir": "differentiation trajectory inference",
    "scVelo": "RNA velocity analysis tool",
    "velocyto": "RNA velocity estimation tool",
    "Scanpy": "single-cell analysis in Python",
    "Seurat": "R toolkit for single-cell genomics",
    "AnnData": "annotated data matrix for single-cell",
    "MuData": "multi-modal annotated data",
    "Squidpy": "spatial single-cell analysis in Python",
    "cell2location": "spatial deconvolution method",
    "RCTD": "robust cell type decomposition",
    "SPOTlight": "spatial transcriptomics deconvolution",
    "BayesSpace": "Bayesian spatial transcriptomics clustering",

    # -- TME / immuno-oncology --
    "TME": "tumor microenvironment",
    "ICB": "immune checkpoint blockade",
    "IO": "immuno-oncology",
    "ICI": "immune checkpoint inhibitor",
    "ADCC": "antibody-dependent cellular cytotoxicity",
    "CDC": "complement-dependent cytotoxicity",
    "ADC": "antibody-drug conjugate",
    "CAR-T": "chimeric antigen receptor T cell",
    "CART": "chimeric antigen receptor T cell",
    "BiTE": "bispecific T-cell engager",
    "TCE": "T-cell engager",
    "ACT": "adoptive cell therapy",
    "TIL therapy": "tumor-infiltrating lymphocyte therapy",
    "CRS": "cytokine release syndrome",
    "ICANS": "immune effector cell-associated neurotoxicity syndrome",
    "MRD": "minimal residual disease",
    "EMT": "epithelial-mesenchymal transition",
    "MET": "mesenchymal-epithelial transition",

    # -- Cancer types --
    "NSCLC": "non-small cell lung cancer",
    "SCLC": "small cell lung cancer",
    "TNBC": "triple-negative breast cancer",
    "CRC": "colorectal cancer",
    "HCC": "hepatocellular carcinoma",
    "RCC": "renal cell carcinoma",
    "GBM": "glioblastoma multiforme",
    "AML": "acute myeloid leukemia",
    "ALL": "acute lymphoblastic leukemia",
    "CLL": "chronic lymphocytic leukemia",
    "CML": "chronic myeloid leukemia",
    "DLBCL": "diffuse large B-cell lymphoma",
    "MM": "multiple myeloma",
    "MDS": "myelodysplastic syndrome",
    "GIST": "gastrointestinal stromal tumor",
    "PDAC": "pancreatic ductal adenocarcinoma",
    "ccRCC": "clear cell renal cell carcinoma",
    "HNSCC": "head and neck squamous cell carcinoma",
    "GEJ": "gastroesophageal junction cancer",

    # -- Foundation models --
    "scGPT": "single-cell generative pre-trained transformer",
    "Geneformer": "Geneformer single-cell foundation model",
    "scFoundation": "large-scale single-cell foundation model",
    "scBERT": "single-cell BERT model",
    "scVI": "single-cell variational inference",
    "scANVI": "single-cell annotation using variational inference",
    "totalVI": "total variational inference for CITE-seq",
    "MultiVI": "multi-omic variational inference",
    "CellTypist": "automated cell type annotation tool",
    "Azimuth": "reference-based single-cell annotation",

    # -- Databases and resources --
    "CellxGene": "CZ CELLxGENE single-cell data portal",
    "HCA": "Human Cell Atlas",
    "GEO": "Gene Expression Omnibus",
    "SRA": "Sequence Read Archive",
    "DepMap": "Cancer Dependency Map",
    "GDSC": "Genomics of Drug Sensitivity in Cancer",
    "CCLE": "Cancer Cell Line Encyclopedia",
    "TCGA": "The Cancer Genome Atlas",

    # -- T cell memory / effector states (new keys) --
    "TRM": "tissue-resident memory T cell",
    "TCM": "central memory T cell",
    "TEM": "effector memory T cell",
    "Teff": "effector T cell",
    "Tex": "exhausted T cell",

    # -- Single-cell technologies (new keys) --
    "scATAC": "single-cell ATAC-seq chromatin accessibility",
    "CyTOF": "cytometry by time of flight mass cytometry",

    # -- Batch correction / integration tools (new keys) --
    "Harmony": "Harmony batch correction integration method",
    "BBKNN": "batch balanced k-nearest neighbors integration",
    "Scanorama": "Scanorama panoramic stitching batch correction",
    "fastMNN": "fast mutual nearest neighbors batch correction",

    # -- Spatial deconvolution methods (new keys) --
    "Tangram": "spatial deconvolution and mapping tool",
    "CIBERSORTx": "bulk RNA-seq deconvolution into cell type fractions",
    "MuSiC": "multi-subject single-cell deconvolution",
    "stereoscope": "spatial transcriptomics deconvolution via scRNA-seq",
    "DestDE": "spatially aware differential expression for deconvolution",

    # -- Sequencing platforms (new keys) --
    "10x Chromium": "10x Genomics Chromium droplet-based single-cell platform",

    # -- Multi-omic assays (new keys) --
    "DOGMA-seq": "simultaneous profiling of chromatin accessibility, RNA, and protein",
    "Multiome": "10x Genomics Multiome joint RNA and ATAC profiling",
    "TotalSeq": "BioLegend TotalSeq antibody-oligo conjugates for CITE-seq",

    # -- Spatial analysis tools (new keys) --
    "SpatialDE": "spatially variable gene detection method",
    "SPARK": "spatial pattern recognition via kernels",
    "ArchR": "analysis of regulatory chromatin in R",
    "Signac": "single-cell chromatin analysis toolkit",
    "chromVAR": "chromatin variability across regions",
    "LIANA": "ligand-receptor analysis framework",
    "LIANA+": "enhanced ligand-receptor interaction analysis",

    # -- Additional cell type states (new keys) --
    "gdT": "gamma-delta T cell",
    "Vd1": "Vdelta1 gamma-delta T cell",
    "Vd2": "Vgamma9Vdelta2 gamma-delta T cell",
    "MAITs": "mucosal-associated invariant T cells",

    # -- Additional database / resource abbreviations (new keys) --
    "PanglaoDB": "PanglaoDB cell type marker database",
    "CellMarker": "CellMarker cell type marker database",
    "CL": "Cell Ontology",
}


# ===================================================================
# SYNONYM MAPS (12 domain-specific synonym maps)
# ===================================================================

CELL_TYPE_MAP: Dict[str, List[str]] = {
    "t_cell": ["T lymphocyte", "T cell", "CD3+ cell", "T-cell", "thymocyte",
               "alpha-beta T cell", "adaptive lymphocyte"],
    "cd8_t": ["cytotoxic T cell", "killer T cell", "CD8+ T cell",
              "CTL", "cytotoxic T lymphocyte", "CD8-positive T cell",
              "effector T cell"],
    "cd4_t": ["helper T cell", "CD4+ T cell", "CD4-positive T cell",
              "T helper cell", "Th cell"],
    "treg": ["regulatory T cell", "T regulatory cell", "suppressor T cell",
             "FOXP3+ T cell", "CD4+CD25+ T cell", "immunoregulatory T cell"],
    "nk_cell": ["natural killer cell", "NK lymphocyte", "innate lymphocyte",
                "CD56+ cell", "cytotoxic innate lymphocyte"],
    "b_cell": ["B lymphocyte", "B cell", "CD19+ cell", "CD20+ cell",
               "immunoglobulin-producing cell precursor"],
    "plasma_cell": ["antibody-secreting cell", "ASC", "plasma B cell",
                    "CD138+ cell", "immunoglobulin-secreting cell",
                    "plasmablast", "long-lived plasma cell"],
    "monocyte": ["CD14+ cell", "myeloid monocyte", "blood monocyte",
                 "classical monocyte", "nonclassical monocyte",
                 "intermediate monocyte"],
    "macrophage": ["tissue macrophage", "histiocyte", "phagocyte",
                   "M1 macrophage", "M2 macrophage", "TAM",
                   "tumor-associated macrophage", "alveolar macrophage",
                   "Kupffer cell"],
    "dendritic_cell": ["DC", "antigen-presenting cell", "APC",
                       "conventional DC", "plasmacytoid DC",
                       "Langerhans cell", "myeloid DC"],
    "fibroblast": ["stromal cell", "mesenchymal cell", "connective tissue cell",
                   "myofibroblast", "cancer-associated fibroblast", "CAF",
                   "activated fibroblast"],
    "endothelial": ["endothelial cell", "vascular cell", "CD31+ cell",
                    "VWF+ cell", "blood vessel cell", "angiogenic cell",
                    "lymphatic endothelial cell"],
    "epithelial": ["epithelial cell", "EpCAM+ cell", "cytokeratin+ cell",
                   "luminal cell", "basal cell", "secretory cell",
                   "ciliated cell", "goblet cell"],
    "neutrophil": ["PMN", "polymorphonuclear cell", "granulocyte",
                   "S100A8+ cell", "N1 neutrophil", "N2 neutrophil",
                   "tumor-associated neutrophil", "TAN"],
    "mast_cell": ["KIT+ granulocyte", "tryptase+ cell", "tissue mast cell",
                  "mucosal mast cell", "connective tissue mast cell"],
    "stem_cell": ["HSC", "hematopoietic stem cell", "CD34+ cell",
                  "progenitor cell", "multipotent progenitor",
                  "long-term HSC", "short-term HSC"],
}

MARKER_MAP: Dict[str, List[str]] = {
    "cd3": ["CD3D", "CD3E", "CD3G", "T-cell marker", "pan-T marker",
            "T-cell receptor complex"],
    "cd8": ["CD8A", "CD8B", "cytotoxic T marker", "MHC-I coreceptor"],
    "cd4": ["CD4 molecule", "helper T marker", "MHC-II coreceptor"],
    "foxp3": ["FOXP3", "forkhead box P3", "Treg marker",
              "regulatory T cell transcription factor"],
    "cd19": ["CD19", "B-cell marker", "B lymphocyte antigen CD19"],
    "cd20": ["MS4A1", "B-cell marker CD20", "rituximab target"],
    "pd1": ["PDCD1", "PD-1", "programmed death 1", "checkpoint receptor",
            "CD279", "immune checkpoint"],
    "pdl1": ["CD274", "PD-L1", "programmed death ligand 1", "B7-H1",
             "immune checkpoint ligand"],
    "checkpoint": ["immune checkpoint", "inhibitory receptor", "PD-1",
                   "CTLA-4", "LAG-3", "TIM-3", "TIGIT"],
    "granzyme": ["GZMA", "GZMB", "GZMK", "cytotoxic granule",
                 "serine protease", "cytotoxicity marker"],
}

SPATIAL_MAP: Dict[str, List[str]] = {
    "visium": ["10x Visium", "spatial transcriptomics", "spot-based",
               "55 micron", "Visium CytAssist", "Visium HD",
               "capture-based spatial"],
    "merfish": ["MERFISH", "Vizgen", "MERSCOPE", "multiplexed FISH",
                "subcellular spatial", "single-molecule FISH"],
    "xenium": ["10x Xenium", "Xenium in situ", "panel-based spatial",
               "in situ hybridization", "Xenium panel"],
    "codex": ["CODEX", "Akoya", "co-detection by indexing",
              "spatial proteomics", "antibody-based spatial",
              "PhenoCycler"],
    "cosmx": ["CosMx", "NanoString", "spatial molecular imager",
              "CosMx SMI", "NanoString spatial"],
    "slide_seq": ["Slide-seq", "Slide-seqV2", "bead-based spatial",
                  "10 micron resolution"],
    "stereo_seq": ["Stereo-seq", "BGI spatial", "spatiotemporal omics",
                   "subcellular spatial"],
    "deconvolution": ["spot deconvolution", "cell2location", "RCTD",
                      "SPOTlight", "Tangram", "BayesSpace",
                      "spatial deconvolution"],
}

TME_MAP: Dict[str, List[str]] = {
    "hot_tumor": ["immune inflamed", "T-cell inflamed", "hot TME",
                  "immune infiltrated", "inflamed tumor",
                  "high TIL density", "PD-L1 positive"],
    "cold_tumor": ["immune desert", "cold TME", "immune excluded",
                   "non-inflamed", "low TIL", "immune cold",
                   "immunologically ignorant"],
    "excluded": ["immune excluded", "stromal barrier", "T cells at margin",
                 "peritumoral immune", "TGF-beta driven exclusion",
                 "fibrotic exclusion"],
    "immunosuppressive": ["immunosuppressive TME", "Treg enriched",
                          "MDSC enriched", "M2 polarized",
                          "immunosuppressive milieu", "IDO high",
                          "arginase high"],
    "immune_checkpoint": ["checkpoint expression", "PD-L1 expression",
                          "checkpoint blockade", "ICB response",
                          "immunotherapy response", "anti-PD-1 response"],
    "tme_profiling": ["tumor microenvironment profiling", "TME classification",
                      "immune contexture", "immune landscape",
                      "immune phenotyping", "tumor immune status"],
}

DRUG_RESPONSE_MAP: Dict[str, List[str]] = {
    "checkpoint_inhibitor": ["anti-PD-1", "anti-PD-L1", "anti-CTLA-4",
                             "pembrolizumab", "nivolumab", "ipilimumab",
                             "atezolizumab", "durvalumab", "avelumab",
                             "immune checkpoint inhibitor", "ICI"],
    "targeted_therapy": ["tyrosine kinase inhibitor", "TKI", "EGFR inhibitor",
                         "BRAF inhibitor", "MEK inhibitor", "CDK4/6 inhibitor",
                         "PARP inhibitor", "mTOR inhibitor"],
    "cell_therapy": ["CAR-T", "chimeric antigen receptor", "adoptive cell therapy",
                     "TIL therapy", "NK cell therapy", "bispecific antibody",
                     "T-cell engager"],
    "drug_sensitivity": ["drug response", "IC50", "dose-response",
                         "pharmacogenomics", "drug resistance",
                         "sensitive cell type", "resistant cell type"],
    "resistance": ["drug resistance", "acquired resistance",
                   "intrinsic resistance", "resistance mechanism",
                   "resistant subpopulation", "resistant clone",
                   "resistance signature"],
}

TRAJECTORY_MAP: Dict[str, List[str]] = {
    "pseudotime": ["pseudotime analysis", "trajectory inference",
                   "cell ordering", "temporal ordering",
                   "differentiation trajectory", "developmental path"],
    "rna_velocity": ["RNA velocity", "spliced unspliced", "scVelo",
                     "velocyto", "dynamical model", "latent time",
                     "transcriptional dynamics"],
    "differentiation": ["cell differentiation", "lineage commitment",
                        "cell fate decision", "progenitor to mature",
                        "developmental trajectory"],
    "exhaustion_trajectory": ["T-cell exhaustion", "exhaustion trajectory",
                              "progenitor exhausted", "terminally exhausted",
                              "TOX expression", "chronic stimulation"],
    "emt": ["epithelial mesenchymal transition", "EMT", "MET",
            "mesenchymal transition", "E-cadherin loss",
            "vimentin expression", "hybrid EMT"],
    "cell_cycle": ["cell cycle phase", "G1 S G2M", "proliferation",
                   "cycling cells", "cell cycle scoring",
                   "S phase", "G2M phase"],
}

TECHNOLOGY_MAP: Dict[str, List[str]] = {
    "droplet": ["droplet-based", "10x Chromium", "Drop-seq", "inDrop",
                "microfluidic", "GEM", "gel bead"],
    "plate_based": ["plate-based", "Smart-seq2", "SMART-seq", "CEL-Seq2",
                    "well-based", "full-length transcript"],
    "combinatorial": ["combinatorial indexing", "sci-RNA-seq", "split-pool",
                      "SPLiT-seq", "sci-RNA-seq3"],
    "multiome": ["multiomics", "multi-omic", "RNA+ATAC", "SHARE-seq",
                 "10x Multiome", "joint profiling", "CITE-seq",
                 "TEA-seq"],
    "perturbation": ["perturbation screen", "Perturb-seq", "CROP-seq",
                     "CRISPR screen", "genetic perturbation",
                     "perturbation atlas"],
    "single_nucleus": ["single-nucleus", "snRNA-seq", "snATAC-seq",
                       "nucleus isolation", "frozen tissue",
                       "nuclear transcriptome"],
}

CANCER_TYPE_MAP: Dict[str, List[str]] = {
    "lung": ["lung cancer", "NSCLC", "SCLC", "lung adenocarcinoma",
             "lung squamous", "pulmonary carcinoma"],
    "breast": ["breast cancer", "TNBC", "triple negative", "ER positive",
               "HER2 positive", "luminal", "basal-like"],
    "melanoma": ["melanoma", "cutaneous melanoma", "uveal melanoma",
                 "BRAF mutant", "NRAS mutant", "skin cancer"],
    "colorectal": ["colorectal cancer", "CRC", "colon cancer",
                   "rectal cancer", "MSI-high", "microsatellite instable"],
    "glioblastoma": ["glioblastoma", "GBM", "brain tumor", "glioma",
                     "high-grade glioma", "IDH-wildtype"],
    "leukemia": ["leukemia", "AML", "ALL", "CLL", "CML",
                 "acute leukemia", "chronic leukemia",
                 "myeloid leukemia", "lymphoblastic leukemia"],
    "lymphoma": ["lymphoma", "DLBCL", "follicular lymphoma",
                 "Hodgkin lymphoma", "non-Hodgkin lymphoma",
                 "mantle cell lymphoma"],
    "myeloma": ["multiple myeloma", "MM", "plasma cell neoplasm",
                "myeloma", "plasmacytoma"],
    "pancreatic": ["pancreatic cancer", "PDAC", "pancreatic ductal",
                   "pancreatic adenocarcinoma"],
    "renal": ["renal cell carcinoma", "RCC", "kidney cancer",
              "clear cell RCC", "ccRCC", "papillary RCC"],
}

TISSUE_MAP: Dict[str, List[str]] = {
    "pbmc": ["PBMC", "peripheral blood", "blood", "circulating cells",
             "whole blood", "buffy coat"],
    "tumor": ["tumor biopsy", "tumor tissue", "neoplasm", "malignancy",
              "tumor sample", "tumor resection", "tumor core"],
    "lymph_node": ["lymph node", "draining lymph node", "sentinel node",
                   "lymphatic tissue", "LN"],
    "bone_marrow": ["bone marrow", "BM aspirate", "marrow",
                    "hematopoietic niche", "BM biopsy"],
    "brain": ["brain tissue", "cerebral cortex", "cerebellum",
              "hippocampus", "CNS tissue", "neural tissue"],
    "lung": ["lung tissue", "airway", "alveolar", "bronchial",
             "pulmonary", "respiratory epithelium"],
    "gut": ["intestine", "colon", "ileum", "jejunum", "duodenum",
            "intestinal mucosa", "GI tract", "gut mucosa"],
    "liver": ["liver tissue", "hepatic", "liver biopsy",
              "liver parenchyma", "hepatic lobule"],
    "skin": ["skin biopsy", "dermis", "epidermis", "cutaneous",
             "skin tissue", "integument"],
    "kidney": ["kidney tissue", "renal", "nephron", "glomerulus",
               "tubular", "renal cortex", "renal medulla"],
}

METHOD_MAP: Dict[str, List[str]] = {
    "clustering": ["cell clustering", "Leiden clustering", "Louvain clustering",
                   "graph-based clustering", "community detection",
                   "Phenograph", "k-means clustering"],
    "annotation": ["cell type annotation", "cell type identification",
                   "cell type classification", "automated annotation",
                   "reference mapping", "label transfer",
                   "CellTypist", "Azimuth", "SingleR"],
    "integration": ["batch integration", "data integration",
                    "batch correction", "Harmony", "scVI",
                    "BBKNN", "Scanorama", "MNN", "fastMNN",
                    "dataset integration"],
    "dimensionality_reduction": ["UMAP", "tSNE", "PCA", "diffusion map",
                                 "force-directed layout", "embedding",
                                 "low-dimensional representation"],
    "differential_expression": ["DE analysis", "differential expression",
                                "DEG analysis", "Wilcoxon test", "MAST",
                                "edgeR", "DESeq2", "pseudobulk DE"],
    "imputation": ["gene imputation", "dropout imputation", "MAGIC",
                   "scImpute", "SAVER", "DCA", "expression recovery",
                   "zero inflation"],
    "doublet_detection": ["doublet removal", "Scrublet", "DoubletFinder",
                          "doublet detection", "multiplet removal",
                          "demultiplexing"],
    "cell_communication": ["cell-cell communication", "ligand-receptor",
                           "CellChat", "CellPhoneDB", "NicheNet",
                           "intercellular signaling", "cell interaction",
                           "LIANA", "cell-cell interaction"],
}

IMMUNE_MAP: Dict[str, List[str]] = {
    "exhaustion": ["T-cell exhaustion", "exhausted T cells", "PD-1 high",
                   "TOX expression", "terminal exhaustion",
                   "progenitor exhausted", "dysfunctional T cell",
                   "chronic activation"],
    "cytotoxicity": ["cytotoxic activity", "granzyme expression",
                     "perforin expression", "cytolytic activity",
                     "killing capacity", "GZMB expression",
                     "effector function"],
    "immune_infiltration": ["immune infiltrate", "immune score",
                            "TIL score", "immune abundance",
                            "immune cell proportion", "ESTIMATE score",
                            "MCP-counter", "xCell"],
    "immune_evasion": ["immune evasion", "immune escape",
                       "antigen loss", "MHC downregulation",
                       "B2M loss", "neoantigen depletion",
                       "immunoediting"],
    "inflammation": ["inflammatory response", "interferon signaling",
                     "IFN-gamma", "TNF signaling", "NF-kB pathway",
                     "pro-inflammatory", "cytokine storm"],
}

GENE_EXPRESSION_MAP: Dict[str, List[str]] = {
    "normalization": ["log normalization", "sctransform", "scran",
                      "Pearson residuals", "library size normalization",
                      "count normalization", "TPM", "CPM"],
    "hvg": ["highly variable genes", "feature selection", "variable genes",
            "top genes", "informative genes", "gene filtering"],
    "scoring": ["gene set scoring", "module scoring", "AddModuleScore",
                "AUCell", "ssGSEA", "pathway scoring",
                "signature scoring", "GSVA"],
    "regulon": ["regulon", "transcription factor activity", "SCENIC",
                "gene regulatory network", "TF target", "GRN inference",
                "regulon activity"],
    "cnv_inference": ["inferred CNV", "InferCNV", "CopyKAT", "Numbat",
                      "copy number from scRNA", "clonal CNV",
                      "aneuploidy inference"],
}

SPATIAL_ANALYSIS_MAP: Dict[str, List[str]] = {
    "spatial_deconvolution": ["cell2location", "Tangram", "RCTD", "stereoscope",
                              "DestDE", "SPOTlight", "BayesSpace",
                              "spatial cell type estimation"],
    "spatial_autocorrelation": ["Moran's I", "Geary's C", "spatial correlation",
                                "spatially variable gene", "SpatialDE",
                                "SPARK", "spatial statistics"],
    "niche_detection": ["neighborhood enrichment", "Ripley's K", "Ripley's L",
                        "spatial niche", "cell neighborhood", "niche clustering",
                        "spatial domain", "tissue compartment",
                        "microenvironment niche"],
}

MULTIOMICS_MAP: Dict[str, List[str]] = {
    "multiome_assay": ["CITE-seq", "Multiome", "SHARE-seq", "DOGMA-seq",
                       "TEA-seq", "ECCITE-seq", "joint RNA+ATAC",
                       "multi-omic single-cell"],
    "chromatin_accessibility": ["ATAC peaks", "open chromatin", "motif enrichment",
                                "transcription factor motif", "chromVAR",
                                "ArchR", "Signac", "peak-to-gene linkage"],
    "protein_surface": ["ADT antibody", "TotalSeq", "antibody-derived tag",
                        "surface protein", "REAP-seq", "CITE-seq protein",
                        "protein quantification", "oligo-conjugated antibody"],
}

# Aggregate all synonym maps for unified lookup
SC_SYNONYMS: Dict[str, Dict[str, List[str]]] = {
    "cell_type": CELL_TYPE_MAP,
    "marker": MARKER_MAP,
    "spatial": SPATIAL_MAP,
    "tme": TME_MAP,
    "drug_response": DRUG_RESPONSE_MAP,
    "trajectory": TRAJECTORY_MAP,
    "technology": TECHNOLOGY_MAP,
    "cancer_type": CANCER_TYPE_MAP,
    "tissue": TISSUE_MAP,
    "method": METHOD_MAP,
    "immune": IMMUNE_MAP,
    "gene_expression": GENE_EXPRESSION_MAP,
    "spatial_analysis": SPATIAL_ANALYSIS_MAP,
    "multiomics": MULTIOMICS_MAP,
}


# ===================================================================
# WORKFLOW TERMS
# ===================================================================

# Maps to SCWorkflowType values defined in src.models.

_WORKFLOW_TERMS: Dict[str, List[str]] = {
    "cell_type_annotation": [
        "cell type", "annotation", "cell identity", "marker gene",
        "cluster identity", "cell classification", "cell label",
        "CellTypist", "Azimuth", "SingleR", "reference mapping",
        "label transfer", "cell ontology", "lineage",
        "CD3", "CD8", "CD4", "FOXP3", "CD19", "CD14", "EPCAM",
        "COL1A1", "PECAM1", "cell atlas",
    ],
    "tme_profiling": [
        "tumor microenvironment", "TME", "immune infiltration",
        "hot tumor", "cold tumor", "immune desert", "immune excluded",
        "immunosuppressive", "TIL", "immune score", "stromal score",
        "immune phenotype", "immune contexture", "checkpoint expression",
        "PD-L1", "immune landscape", "immune composition",
        "T cell infiltration", "macrophage polarization",
    ],
    "drug_response": [
        "drug response", "drug sensitivity", "drug resistance",
        "IC50", "pharmacogenomics", "treatment response",
        "checkpoint inhibitor", "targeted therapy", "immunotherapy",
        "resistant clone", "sensitive cell type", "combination therapy",
        "synergy", "GDSC", "DepMap", "therapeutic vulnerability",
    ],
    "subclonal_architecture": [
        "subclone", "clonal", "clone", "subpopulation",
        "clonal evolution", "tumor heterogeneity", "InferCNV",
        "CopyKAT", "copy number", "CNV", "clonal architecture",
        "phylogeny", "driver mutation", "clonal fitness",
        "tumor evolution", "intratumoral heterogeneity",
    ],
    "spatial_niche": [
        "spatial", "niche", "spatial transcriptomics", "Visium",
        "MERFISH", "Xenium", "CODEX", "CosMx", "tissue architecture",
        "spatial pattern", "cell neighborhood", "spatial domain",
        "spatial cluster", "cell-cell proximity", "tissue region",
        "spatial deconvolution", "spatially variable gene",
    ],
    "trajectory_analysis": [
        "trajectory", "pseudotime", "differentiation", "lineage",
        "RNA velocity", "cell fate", "branching point",
        "developmental path", "Monocle", "scVelo", "PAGA",
        "Palantir", "diffusion pseudotime", "cell ordering",
        "transition", "progenitor", "terminal state",
    ],
    "ligand_receptor": [
        "ligand-receptor", "cell communication", "cell interaction",
        "CellChat", "CellPhoneDB", "NicheNet", "signaling",
        "receptor", "ligand", "intercellular", "crosstalk",
        "paracrine", "juxtacrine", "cell-cell signaling",
        "communication network", "interaction score",
    ],
    "biomarker_discovery": [
        "biomarker", "diagnostic marker", "prognostic marker",
        "predictive biomarker", "surface marker", "drug target",
        "differential expression", "fold change", "specificity",
        "sensitivity", "AUROC", "clinical assay",
        "companion diagnostic", "pharmacodynamic marker",
    ],
    "cart_target_validation": [
        "CAR-T", "chimeric antigen receptor", "target validation",
        "antigen expression", "on-target off-tumor", "toxicity",
        "antigen loss", "escape variant", "CD19 CAR-T", "BCMA",
        "target antigen", "surface protein", "tumor antigen",
        "co-expression", "dual targeting", "safety profile",
    ],
    "treatment_monitoring": [
        "treatment monitoring", "longitudinal", "treatment response",
        "pre-treatment", "post-treatment", "on-treatment",
        "baseline", "timepoint", "compositional shift",
        "emerging clone", "resistance emergence", "MRD",
        "minimal residual disease", "immune dynamics",
        "clonal evolution under treatment",
    ],
    "general": [
        "single-cell", "single cell", "scRNA-seq", "single-nucleus",
        "transcriptomics", "gene expression", "cell population",
        "clustering", "UMAP", "embedding", "quality control",
    ],
}


# ===================================================================
# QUERY EXPANDER CLASS
# ===================================================================


class QueryExpander:
    """Expand single-cell queries with aliases, synonyms, and workflow terms.

    The expander performs three operations:
    1. **Entity detection** -- identifies known abbreviations, technology
       names, and cell type aliases in the query text and resolves them
       to canonical forms.
    2. **Synonym expansion** -- augments the query with domain-specific
       synonyms drawn from the 12 synonym maps.
    3. **Workflow term injection** -- adds terms relevant to the detected
       or specified SCWorkflowType to improve recall.
    """

    def __init__(
        self,
        aliases: Optional[Dict[str, str]] = None,
        synonyms: Optional[Dict[str, Dict[str, List[str]]]] = None,
        workflow_terms: Optional[Dict[str, List[str]]] = None,
        max_expansion_terms: int = 30,
    ) -> None:
        self._aliases = aliases or ENTITY_ALIASES
        self._synonyms = synonyms or SC_SYNONYMS
        self._workflow_terms = workflow_terms or _WORKFLOW_TERMS
        self._max_expansion_terms = max_expansion_terms

        # Pre-compute a case-insensitive lookup for aliases
        self._alias_lower: Dict[str, str] = {
            k.lower(): v for k, v in self._aliases.items()
        }

        # Build a flat reverse index: synonym term -> category
        self._synonym_index: Dict[str, str] = {}
        for domain, category_map in self._synonyms.items():
            for category, terms in category_map.items():
                for term in terms:
                    self._synonym_index[term.lower()] = (
                        f"{domain}.{category}"
                    )

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def expand(
        self, query: str, workflow: Optional[str] = None
    ) -> List[str]:
        """Return a list of expansion terms for the given query.

        Parameters
        ----------
        query : str
            Free-text single-cell analysis query.
        workflow : str, optional
            SCWorkflowType value (e.g. ``"cell_type_annotation"``).
            If ``None``, the expander attempts to infer the workflow
            from query content.

        Returns
        -------
        list[str]
            Deduplicated list of expansion terms, capped at
            ``max_expansion_terms``.
        """
        terms: list[str] = []
        query_lower = query.lower()

        # 1. Resolve entities / aliases
        detected = self.detect_entities(query)
        for _alias, canonical in detected.items():
            terms.append(canonical)

        # 2. Synonym expansion
        for domain, category_map in self._synonyms.items():
            for category, syns in category_map.items():
                # Check if query mentions category key or any synonym
                if category.lower() in query_lower:
                    terms.extend(syns[:5])  # top 5 per match
                    continue
                for syn in syns:
                    if syn.lower() in query_lower:
                        terms.extend(syns[:5])
                        break

        # 3. Workflow-specific terms
        wf_key = workflow or self._infer_workflow(query_lower)
        wf_terms = self.get_workflow_terms(wf_key)
        terms.extend(wf_terms)

        # Deduplicate while preserving order, cap at limit
        seen: set[str] = set()
        unique: list[str] = []
        for t in terms:
            t_lower = t.lower()
            if t_lower not in seen and t_lower not in query_lower:
                seen.add(t_lower)
                unique.append(t)
        return unique[: self._max_expansion_terms]

    def detect_entities(self, query: str) -> Dict[str, str]:
        """Detect abbreviations and technology names in the query.

        Returns
        -------
        dict[str, str]
            Mapping of matched alias -> canonical expansion.
        """
        detected: Dict[str, str] = {}
        _STOP_WORDS = {
            "a", "an", "and", "are", "as", "at", "be", "by", "for",
            "from", "had", "has", "have", "he", "her", "his", "how",
            "i", "if", "in", "is", "it", "its", "my", "no", "not",
            "of", "on", "or", "our", "she", "so", "than", "that",
            "the", "then", "they", "this", "to", "up", "us", "was",
            "we", "what", "when", "who", "will", "with", "you",
        }
        # Tokenize into words and multi-word spans
        tokens = re.findall(r"[A-Za-z0-9\-'+]+", query)
        for token in tokens:
            token_lower = token.lower()
            if token_lower in _STOP_WORDS:
                continue
            # Exact-case match first (case-sensitive abbreviations
            # like MERFISH, scRNA-seq, PD-1)
            if token in self._aliases:
                detected[token] = self._aliases[token]
            elif token_lower in self._alias_lower:
                detected[token] = self._alias_lower[token_lower]
        return detected

    def get_workflow_terms(
        self, workflow: Optional[str] = None
    ) -> List[str]:
        """Return search terms associated with the given workflow.

        Parameters
        ----------
        workflow : str, optional
            SCWorkflowType value.  Falls back to ``"general"``
            if not recognized.

        Returns
        -------
        list[str]
            Search terms for the workflow.
        """
        if workflow is None:
            return self._workflow_terms.get("general", [])
        return self._workflow_terms.get(
            workflow, self._workflow_terms.get("general", [])
        )

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    def _infer_workflow(self, query_lower: str) -> str:
        """Best-effort workflow inference from query text.

        Scores each workflow by counting how many of its terms appear
        in the query and picks the highest-scoring one.
        """
        best_wf = "general"
        best_score = 0
        for wf_key, terms in self._workflow_terms.items():
            if wf_key == "general":
                continue
            score = sum(1 for t in terms if t.lower() in query_lower)
            if score > best_score:
                best_score = score
                best_wf = wf_key
        return best_wf if best_score > 0 else "general"
