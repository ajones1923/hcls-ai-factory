"""Single-Cell Intelligence Agent -- Domain Knowledge Base.

Comprehensive single-cell genomics knowledge covering 30+ cell types,
tumor microenvironment profiles, drug sensitivity databases, spatial
transcriptomics platforms, immune signatures, foundation models,
GPU benchmarks, and canonical marker gene databases for single-cell
RAG-based analysis and clinical decision support.

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict


# ===================================================================
# KNOWLEDGE BASE VERSION
# ===================================================================

KNOWLEDGE_VERSION: Dict[str, Any] = {
    "version": "3.0.0",
    "last_updated": "2026-03-22",
    "revision_notes": "Major expansion -- 44 cell types with Cell Ontology IDs, "
                      "4 TME profiles with clinical trials, 30 drugs with IC50 "
                      "ranges, 4 spatial platforms, 75 marker genes, 10 immune "
                      "signatures, 3 foundation models, GPU benchmarks, "
                      "25 ligand-receptor pairs, 12 cancer TME atlas profiles.",
    "sources": [
        "Human Cell Atlas",
        "CellMarker 2.0",
        "Cell Ontology (CL)",
        "Tabula Sapiens",
        "Tabula Muris",
        "PanglaoDB",
        "CancerSEA",
        "GDSC (Genomics of Drug Sensitivity in Cancer)",
        "DepMap",
        "10x Genomics Spatial Protocols",
        "Vizgen MERFISH Documentation",
        "scGPT / Geneformer / scFoundation Publications",
        "PubMed / MEDLINE",
        "TISCH2 (Tumor Immune Single-Cell Hub)",
        "CellTalkDB",
        "CellPhoneDB",
        "NicheNet",
        "CIBERSORTx",
        "TCGA Pan-Cancer Atlas",
    ],
    "counts": {
        "cell_types": 44,
        "tme_profiles": 4,
        "drugs": 30,
        "spatial_platforms": 4,
        "marker_genes": 75,
        "immune_signatures": 10,
        "foundation_models": 3,
        "gpu_benchmarks": 4,
        "ligand_receptor_pairs": 25,
        "cancer_tme_profiles": 12,
    },
}


# ===================================================================
# CELL TYPE ATLAS (32 cell types)
# ===================================================================

CELL_TYPE_ATLAS: Dict[str, Dict[str, Any]] = {
    "T_cell": {
        "markers": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
        "cell_ontology_id": "CL:0000084",
        "description": "Mature T lymphocyte expressing the T-cell receptor complex. "
                       "Central to adaptive immunity including cytotoxic killing, "
                       "helper functions, and immune regulation.",
        "tissues": ["blood", "lymph_node", "spleen", "thymus", "bone_marrow",
                    "tumor", "gut", "lung", "skin"],
        "subtypes": ["CD4_T", "CD8_T", "Treg", "gamma_delta_T", "NKT", "MAIT"],
    },
    "CD8_T": {
        "markers": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG"],
        "cell_ontology_id": "CL:0000625",
        "description": "Cytotoxic T lymphocyte that recognizes MHC class I-restricted "
                       "antigens and kills target cells via perforin/granzyme pathway.",
        "tissues": ["blood", "lymph_node", "spleen", "tumor", "lung", "liver"],
        "subtypes": ["naive_CD8", "effector_CD8", "memory_CD8", "exhausted_CD8",
                     "tissue_resident_memory_CD8"],
    },
    "CD4_T": {
        "markers": ["CD4", "IL7R", "TCF7", "LEF1", "CCR7"],
        "cell_ontology_id": "CL:0000624",
        "description": "Helper T lymphocyte that coordinates adaptive immune responses "
                       "through cytokine secretion and costimulatory interactions.",
        "tissues": ["blood", "lymph_node", "spleen", "thymus", "gut", "tumor"],
        "subtypes": ["Th1", "Th2", "Th17", "Tfh", "naive_CD4", "memory_CD4"],
    },
    "Treg": {
        "markers": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
        "cell_ontology_id": "CL:0000815",
        "description": "Regulatory T cell that suppresses immune responses and "
                       "maintains self-tolerance. High expression of FOXP3 and CD25 "
                       "(IL2RA). Critical in autoimmunity and tumor immune evasion.",
        "tissues": ["blood", "lymph_node", "spleen", "tumor", "gut", "skin",
                    "adipose_tissue"],
        "subtypes": ["natural_Treg", "induced_Treg", "effector_Treg",
                     "tissue_resident_Treg"],
    },
    "gamma_delta_T": {
        "markers": ["TRGC1", "TRGC2", "TRDC", "TRDV2", "NKG7"],
        "cell_ontology_id": "CL:0000798",
        "description": "Unconventional T cell expressing gamma-delta TCR with innate-like "
                       "recognition of stress antigens and MHC-unrestricted cytotoxicity.",
        "tissues": ["blood", "skin", "gut", "lung", "liver"],
        "subtypes": ["Vdelta1", "Vdelta2", "Vdelta3"],
    },
    "B_cell": {
        "markers": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"],
        "cell_ontology_id": "CL:0000236",
        "description": "B lymphocyte responsible for humoral immunity through "
                       "antibody production. Key component of adaptive immune system.",
        "tissues": ["blood", "lymph_node", "spleen", "bone_marrow", "tumor",
                    "tonsil", "gut"],
        "subtypes": ["naive_B", "memory_B", "germinal_center_B", "marginal_zone_B",
                     "transitional_B"],
    },
    "Plasma": {
        "markers": ["SDC1", "MZB1", "XBP1", "JCHAIN", "IGHG1"],
        "cell_ontology_id": "CL:0000786",
        "description": "Terminally differentiated B cell specialized for high-rate "
                       "immunoglobulin secretion. Identified by CD138 (SDC1) expression.",
        "tissues": ["bone_marrow", "lymph_node", "spleen", "gut", "tumor"],
        "subtypes": ["short_lived_plasma", "long_lived_plasma", "plasmablast"],
    },
    "NK": {
        "markers": ["KLRD1", "NKG7", "NCAM1", "KLRF1", "GNLY"],
        "cell_ontology_id": "CL:0000623",
        "description": "Natural killer cell providing innate cytotoxic surveillance. "
                       "Kills virus-infected and tumor cells without prior sensitization.",
        "tissues": ["blood", "spleen", "liver", "lung", "tumor", "uterus"],
        "subtypes": ["CD56bright_NK", "CD56dim_NK", "adaptive_NK",
                     "tissue_resident_NK"],
    },
    "Monocyte": {
        "markers": ["CD14", "LYZ", "S100A9", "VCAN", "FCN1"],
        "cell_ontology_id": "CL:0000576",
        "description": "Myeloid cell circulating in blood that differentiates into "
                       "macrophages or dendritic cells upon tissue entry.",
        "tissues": ["blood", "bone_marrow", "spleen"],
        "subtypes": ["classical_CD14", "intermediate_CD14_CD16",
                     "nonclassical_CD16"],
    },
    "Macrophage": {
        "markers": ["CD68", "MARCO", "CSF1R", "MRC1", "MSR1"],
        "cell_ontology_id": "CL:0000235",
        "description": "Tissue-resident phagocyte derived from monocytes or yolk-sac "
                       "progenitors. Performs phagocytosis, antigen presentation, and "
                       "tissue remodeling.",
        "tissues": ["lung", "liver", "spleen", "brain", "bone", "gut",
                    "adipose_tissue", "tumor"],
        "subtypes": ["M1_proinflammatory", "M2_anti_inflammatory",
                     "tumor_associated_macrophage", "alveolar_macrophage",
                     "Kupffer_cell", "microglia"],
    },
    "DC": {
        "markers": ["CLEC9A", "CD1C", "FCER1A", "BATF3", "IRF8"],
        "cell_ontology_id": "CL:0000451",
        "description": "Professional antigen-presenting cell bridging innate and "
                       "adaptive immunity. Captures, processes, and presents antigens "
                       "to T cells via MHC molecules.",
        "tissues": ["blood", "lymph_node", "skin", "lung", "gut", "tumor"],
        "subtypes": ["cDC1", "cDC2", "pDC", "moDC", "Langerhans_cell"],
    },
    "pDC": {
        "markers": ["CLEC4C", "IL3RA", "TCF4", "IRF7", "LILRA4"],
        "cell_ontology_id": "CL:0000784",
        "description": "Plasmacytoid dendritic cell specialized for type I interferon "
                       "production in response to viral nucleic acids.",
        "tissues": ["blood", "lymph_node", "bone_marrow", "tumor"],
        "subtypes": ["activated_pDC", "tolerogenic_pDC"],
    },
    "Neutrophil": {
        "markers": ["CSF3R", "S100A8", "S100A12", "CXCR2", "FCGR3B"],
        "cell_ontology_id": "CL:0000775",
        "description": "Most abundant granulocyte in blood providing first-line "
                       "defense via phagocytosis, degranulation, and NETosis.",
        "tissues": ["blood", "bone_marrow", "spleen", "tumor", "lung", "gut"],
        "subtypes": ["N1_antitumor", "N2_protumor", "low_density_neutrophil"],
    },
    "Mast_cell": {
        "markers": ["KIT", "CPA3", "TPSAB1", "TPSB2", "HDC"],
        "cell_ontology_id": "CL:0000097",
        "description": "Tissue-resident granulocyte containing histamine and "
                       "heparin granules. Central to allergic responses and "
                       "tissue homeostasis.",
        "tissues": ["skin", "gut", "lung", "connective_tissue", "tumor"],
        "subtypes": ["mucosal_mast_cell", "connective_tissue_mast_cell"],
    },
    "Basophil": {
        "markers": ["CLC", "GATA2", "HDC", "IL4", "CCR3"],
        "cell_ontology_id": "CL:0000767",
        "description": "Rare circulating granulocyte involved in Th2 immunity "
                       "and parasitic defense via IL-4 and IL-13 production.",
        "tissues": ["blood", "bone_marrow"],
        "subtypes": [],
    },
    "Eosinophil": {
        "markers": ["CCR3", "SIGLEC8", "EPX", "PRG2", "CLC"],
        "cell_ontology_id": "CL:0000771",
        "description": "Granulocyte involved in allergic inflammation and "
                       "anti-parasitic defense. Contains cytotoxic granule proteins.",
        "tissues": ["blood", "bone_marrow", "gut", "lung", "skin"],
        "subtypes": ["tissue_eosinophil", "inflammatory_eosinophil"],
    },
    "Fibroblast": {
        "markers": ["COL1A1", "DCN", "LUM", "VIM", "PDGFRA"],
        "cell_ontology_id": "CL:0000057",
        "description": "Mesenchymal stromal cell producing extracellular matrix "
                       "components. Key in wound healing, fibrosis, and tissue "
                       "architecture.",
        "tissues": ["skin", "lung", "heart", "liver", "gut", "tumor",
                    "connective_tissue"],
        "subtypes": ["myofibroblast", "cancer_associated_fibroblast",
                     "adventitial_fibroblast", "reticular_fibroblast"],
    },
    "Endothelial": {
        "markers": ["PECAM1", "VWF", "CDH5", "ERG", "KDR"],
        "cell_ontology_id": "CL:0000115",
        "description": "Cell lining blood and lymphatic vessels. Regulates vascular "
                       "permeability, angiogenesis, and immune cell trafficking.",
        "tissues": ["blood_vessel", "heart", "lung", "liver", "brain", "kidney",
                    "tumor"],
        "subtypes": ["arterial_endothelial", "venous_endothelial",
                     "capillary_endothelial", "lymphatic_endothelial",
                     "tip_cell", "stalk_cell"],
    },
    "Epithelial": {
        "markers": ["EPCAM", "KRT18", "KRT19", "CDH1", "KRT8"],
        "cell_ontology_id": "CL:0000066",
        "description": "Cell forming epithelial barriers in skin, gut, and organs. "
                       "Provides structural integrity and selective permeability.",
        "tissues": ["skin", "gut", "lung", "kidney", "liver", "breast",
                    "prostate", "pancreas"],
        "subtypes": ["basal_epithelial", "luminal_epithelial",
                     "ciliated_epithelial", "secretory_epithelial",
                     "club_cell", "goblet_cell",
                     "AT1_pneumocyte", "AT2_pneumocyte"],
    },
    "Hepatocyte": {
        "markers": ["ALB", "APOB", "CYP3A4", "HNF4A", "SERPINA1"],
        "cell_ontology_id": "CL:0000182",
        "description": "Parenchymal liver cell performing metabolic, synthetic, "
                       "and detoxification functions. Produces albumin, clotting "
                       "factors, and bile components.",
        "tissues": ["liver"],
        "subtypes": ["periportal_hepatocyte", "pericentral_hepatocyte",
                     "midlobular_hepatocyte"],
    },
    "Neuron": {
        "markers": ["MAP2", "RBFOX3", "SYP", "SNAP25", "NEFL"],
        "cell_ontology_id": "CL:0000540",
        "description": "Electrically excitable cell transmitting information via "
                       "synaptic connections. Forms the functional unit of the "
                       "nervous system.",
        "tissues": ["brain", "spinal_cord", "peripheral_nerve", "retina",
                    "enteric_nervous_system"],
        "subtypes": ["excitatory_neuron", "inhibitory_neuron",
                     "dopaminergic_neuron", "motor_neuron", "sensory_neuron",
                     "interneuron", "pyramidal_neuron", "Purkinje_cell"],
    },
    "Astrocyte": {
        "markers": ["GFAP", "AQP4", "S100B", "SLC1A3", "ALDH1L1"],
        "cell_ontology_id": "CL:0000127",
        "description": "Star-shaped glial cell supporting neuronal metabolism, "
                       "blood-brain barrier integrity, and synaptic function.",
        "tissues": ["brain", "spinal_cord", "retina"],
        "subtypes": ["protoplasmic_astrocyte", "fibrous_astrocyte",
                     "reactive_astrocyte"],
    },
    "Oligodendrocyte": {
        "markers": ["MBP", "PLP1", "MOG", "MAG", "OLIG2"],
        "cell_ontology_id": "CL:0000128",
        "description": "Glial cell forming myelin sheaths around CNS axons. "
                       "Essential for saltatory conduction and axonal support.",
        "tissues": ["brain", "spinal_cord"],
        "subtypes": ["mature_oligodendrocyte",
                     "oligodendrocyte_precursor_cell"],
    },
    "Microglia": {
        "markers": ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "ITGAM"],
        "cell_ontology_id": "CL:0000129",
        "description": "Brain-resident macrophage originating from yolk-sac "
                       "progenitors. Surveys and maintains CNS homeostasis, "
                       "performs synaptic pruning.",
        "tissues": ["brain", "spinal_cord", "retina"],
        "subtypes": ["homeostatic_microglia", "disease_associated_microglia",
                     "activated_microglia"],
    },
    "Cardiomyocyte": {
        "markers": ["TNNT2", "MYH6", "MYH7", "ACTC1", "MYL2"],
        "cell_ontology_id": "CL:0000746",
        "description": "Contractile muscle cell of the heart responsible for "
                       "rhythmic contraction. Connected via intercalated discs.",
        "tissues": ["heart"],
        "subtypes": ["atrial_cardiomyocyte", "ventricular_cardiomyocyte",
                     "pacemaker_cell", "Purkinje_fiber_cell"],
    },
    "Smooth_muscle": {
        "markers": ["ACTA2", "MYH11", "TAGLN", "CNN1", "DES"],
        "cell_ontology_id": "CL:0000192",
        "description": "Non-striated muscle cell found in vessel walls, airways, "
                       "and hollow organs. Mediates involuntary contraction.",
        "tissues": ["blood_vessel", "gut", "bladder", "uterus", "airway"],
        "subtypes": ["vascular_smooth_muscle", "visceral_smooth_muscle",
                     "airway_smooth_muscle"],
    },
    "Adipocyte": {
        "markers": ["ADIPOQ", "FABP4", "LEP", "PPARG", "PLIN1"],
        "cell_ontology_id": "CL:0000136",
        "description": "Lipid-storing cell with endocrine function producing "
                       "adipokines. Central to energy metabolism and thermogenesis.",
        "tissues": ["adipose_tissue", "bone_marrow", "breast"],
        "subtypes": ["white_adipocyte", "brown_adipocyte", "beige_adipocyte"],
    },
    "Pericyte": {
        "markers": ["PDGFRB", "CSPG4", "RGS5", "NOTCH3", "ACTA2"],
        "cell_ontology_id": "CL:0000669",
        "description": "Mural cell wrapped around capillary endothelium regulating "
                       "blood flow, vessel stability, and blood-brain barrier "
                       "integrity.",
        "tissues": ["brain", "lung", "kidney", "retina", "muscle"],
        "subtypes": ["capillary_pericyte", "arteriolar_pericyte",
                     "venular_pericyte"],
    },
    "Mesenchymal_stem_cell": {
        "markers": ["ENG", "THY1", "NT5E", "ITGB1", "ALCAM"],
        "cell_ontology_id": "CL:0000134",
        "description": "Multipotent stromal cell with capacity to differentiate "
                       "into osteoblasts, chondrocytes, and adipocytes.",
        "tissues": ["bone_marrow", "adipose_tissue", "umbilical_cord",
                    "dental_pulp"],
        "subtypes": ["bone_marrow_MSC", "adipose_derived_MSC",
                     "umbilical_cord_MSC"],
    },
    "HSC": {
        "markers": ["CD34", "KIT", "THY1", "PROM1", "CRHBP"],
        "cell_ontology_id": "CL:0000037",
        "description": "Hematopoietic stem cell residing in bone marrow with "
                       "self-renewal capacity and ability to generate all blood "
                       "cell lineages.",
        "tissues": ["bone_marrow", "fetal_liver", "umbilical_cord_blood"],
        "subtypes": ["long_term_HSC", "short_term_HSC",
                     "multipotent_progenitor"],
    },
    "Erythrocyte_precursor": {
        "markers": ["GYPA", "HBB", "HBA1", "KLF1", "TFRC"],
        "cell_ontology_id": "CL:0000764",
        "description": "Erythroid lineage precursor undergoing hemoglobin synthesis "
                       "and enucleation to produce mature red blood cells.",
        "tissues": ["bone_marrow", "fetal_liver"],
        "subtypes": ["proerythroblast", "basophilic_erythroblast",
                     "polychromatic_erythroblast",
                     "orthochromatic_erythroblast"],
    },
    "Megakaryocyte": {
        "markers": ["ITGA2B", "GP9", "PF4", "PPBP", "TUBB1"],
        "cell_ontology_id": "CL:0000556",
        "description": "Large polyploid bone marrow cell that produces platelets "
                       "via cytoplasmic fragmentation.",
        "tissues": ["bone_marrow", "lung"],
        "subtypes": ["immature_megakaryocyte", "mature_megakaryocyte"],
    },
    "Plasmacytoid_DC": {
        "markers": ["CLEC4C", "IL3RA", "IRF7", "TCF4", "LILRA4"],
        "cell_ontology_id": "CL:0000784",
        "description": "Type I interferon producers specialized for antiviral "
                       "defense. Sense viral nucleic acids via TLR7/TLR9 and "
                       "produce large quantities of IFN-alpha/beta.",
        "tissues": ["blood", "lymph_node", "bone_marrow", "tumor"],
        "subtypes": ["activated_pDC", "tolerogenic_pDC"],
    },
    "cDC1": {
        "markers": ["CLEC9A", "XCR1", "BATF3", "IRF8", "THBD"],
        "cell_ontology_id": "CL:0002394",
        "description": "Cross-presentation specialist dendritic cell that excels "
                       "at presenting exogenous antigens on MHC class I to prime "
                       "CD8+ T-cell responses. Critical for anti-tumor immunity.",
        "tissues": ["lymph_node", "spleen", "lung", "tumor", "skin"],
        "subtypes": ["migratory_cDC1", "resident_cDC1"],
    },
    "cDC2": {
        "markers": ["CD1C", "FCER1A", "CLEC10A", "IRF4", "ITGAX"],
        "cell_ontology_id": "CL:0002399",
        "description": "Conventional dendritic cell subtype specializing in CD4 "
                       "T-cell activation through MHC class II antigen "
                       "presentation. Promotes Th1, Th2, and Th17 polarization.",
        "tissues": ["blood", "lymph_node", "skin", "lung", "gut", "tumor"],
        "subtypes": ["cDC2A", "cDC2B"],
    },
    "Innate_lymphoid_cell_1": {
        "markers": ["TBX21", "IFNG", "IL12RB2", "NCR1", "EOMES"],
        "cell_ontology_id": "CL:0001069",
        "description": "Tissue-resident innate lymphoid cell producing IFN-gamma "
                       "for defense against intracellular pathogens. Innate "
                       "counterpart of Th1 cells.",
        "tissues": ["gut", "liver", "uterus", "salivary_gland"],
        "subtypes": ["ILC1", "ieILC1"],
    },
    "Innate_lymphoid_cell_2": {
        "markers": ["GATA3", "IL13", "IL5", "IL33R", "KLRG1"],
        "cell_ontology_id": "CL:0001070",
        "description": "Tissue-resident innate lymphoid cell driving type 2 "
                       "inflammation through IL-5 and IL-13 production. Key in "
                       "allergic responses and anti-helminth immunity.",
        "tissues": ["lung", "gut", "skin", "adipose_tissue"],
        "subtypes": ["natural_ILC2", "inflammatory_ILC2"],
    },
    "Innate_lymphoid_cell_3": {
        "markers": ["RORC", "IL17A", "IL22", "NCR1", "AHR"],
        "cell_ontology_id": "CL:0001071",
        "description": "Tissue-resident innate lymphoid cell supporting mucosal "
                       "immunity through IL-17 and IL-22 production. Maintains "
                       "epithelial barrier integrity.",
        "tissues": ["gut", "lung", "skin", "tonsil"],
        "subtypes": ["NCR_positive_ILC3", "NCR_negative_ILC3", "LTi_cell"],
    },
    "Gamma_delta_T": {
        "markers": ["TRGV9", "TRDV2", "NKG7", "TRDC", "GNLY"],
        "cell_ontology_id": "CL:0000798",
        "description": "Unconventional T cell expressing gamma-delta TCR that "
                       "bridges innate and adaptive immunity. Recognizes "
                       "phosphoantigens and stress ligands without MHC "
                       "restriction.",
        "tissues": ["blood", "skin", "gut", "lung", "liver", "tumor"],
        "subtypes": ["Vgamma9Vdelta2", "Vdelta1", "Vdelta3"],
    },
    "MAIT_cell": {
        "markers": ["TRAV1-2", "SLC4A10", "KLRB1", "IL18R1", "ZBTB16"],
        "cell_ontology_id": "CL:0000940",
        "description": "Mucosal-associated invariant T cell recognizing microbial "
                       "riboflavin metabolites presented by MR1. Abundant in "
                       "liver and mucosal tissues with rapid effector function.",
        "tissues": ["blood", "liver", "gut", "lung"],
        "subtypes": ["MAIT1", "MAIT17"],
    },
    "Erythroid_progenitor": {
        "markers": ["GYPA", "KLF1", "HBB", "HBA1", "TFRC"],
        "cell_ontology_id": "CL:0000764",
        "description": "Red blood cell development progenitor undergoing "
                       "hemoglobin synthesis and progressive nuclear condensation "
                       "en route to mature erythrocyte formation.",
        "tissues": ["bone_marrow", "fetal_liver"],
        "subtypes": ["BFU-E", "CFU-E", "proerythroblast", "erythroblast"],
    },
    "Schwann_cell": {
        "markers": ["MPZ", "PMP22", "SOX10", "MBP", "EGR2"],
        "cell_ontology_id": "CL:0002573",
        "description": "Peripheral nerve support cell forming myelin sheaths "
                       "around peripheral axons. Essential for saltatory "
                       "conduction and nerve regeneration.",
        "tissues": ["peripheral_nerve", "skin", "gut"],
        "subtypes": ["myelinating_Schwann", "non_myelinating_Schwann",
                     "repair_Schwann"],
    },
    "Podocyte": {
        "markers": ["NPHS1", "NPHS2", "WT1", "SYNPO", "PODXL"],
        "cell_ontology_id": "CL:0000653",
        "description": "Highly specialized kidney glomerular epithelial cell with "
                       "foot processes forming the slit diaphragm for selective "
                       "blood filtration. Loss causes proteinuria.",
        "tissues": ["kidney"],
        "subtypes": [],
    },
    "Goblet_cell": {
        "markers": ["MUC2", "TFF3", "SPDEF", "FCGBP", "CLCA1"],
        "cell_ontology_id": "CL:0000160",
        "description": "Mucus-secreting epithelial cell lining the respiratory "
                       "and intestinal tracts. Produces gel-forming mucins that "
                       "form a protective barrier against pathogens.",
        "tissues": ["gut", "lung", "conjunctiva"],
        "subtypes": ["intestinal_goblet", "airway_goblet"],
    },
}


# ===================================================================
# TUMOR MICROENVIRONMENT (TME) PROFILES
# ===================================================================

TME_PROFILES: Dict[str, Dict[str, Any]] = {
    "hot": {
        "description": "Immune-inflamed tumor with high T-cell infiltration "
                       "throughout the tumor parenchyma. Associated with better "
                       "checkpoint inhibitor response.",
        "signature_genes": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "CXCL9",
                           "CXCL10", "CXCL11", "STAT1", "IDO1", "GBP1",
                           "HLA-DRA"],
        "immune_infiltration_thresholds": {
            "CD8_T_fraction": ">0.15",
            "cytotoxic_score": ">0.6",
            "IFNg_signature": ">0.5",
            "TIL_density": ">500_per_mm2",
        },
        "treatment_implications": [
            "High likelihood of checkpoint inhibitor response",
            "Anti-PD-1/PD-L1 monotherapy often effective",
            "Consider combination with anti-CTLA-4 for further benefit",
            "Monitor for immune-related adverse events",
        ],
        "clinical_trials": [
            "KEYNOTE-001", "KEYNOTE-024", "CheckMate 067",
            "IMvigor210", "POPLAR",
        ],
    },
    "cold": {
        "description": "Immune-desert tumor with minimal immune cell infiltration. "
                       "Requires strategies to prime anti-tumor immunity.",
        "signature_genes": ["CTNNB1", "WNT2", "DKK1", "PTEN_loss_signature",
                           "MYC", "CDK4", "CCND1"],
        "immune_infiltration_thresholds": {
            "CD8_T_fraction": "<0.02",
            "cytotoxic_score": "<0.1",
            "IFNg_signature": "<0.1",
            "TIL_density": "<50_per_mm2",
        },
        "treatment_implications": [
            "Checkpoint inhibitor monotherapy unlikely to work",
            "Consider oncolytic virus or radiation to prime immunity",
            "Adoptive cell therapy (CAR-T, TIL therapy) may bypass",
            "Combination with targeted therapy or chemotherapy",
            "Cancer vaccine strategies to generate neoantigen response",
        ],
        "clinical_trials": [
            "MASTERKEY-265 (T-VEC + pembro)", "LEAP trials",
            "RELATIVITY-047",
        ],
    },
    "excluded": {
        "description": "Immune-excluded tumor where immune cells are restricted "
                       "to the tumor periphery/stroma and do not penetrate the "
                       "tumor parenchyma.",
        "signature_genes": ["TGFB1", "TGFB2", "COL1A1", "FAP", "ACTA2",
                           "FN1", "POSTN", "VEGFA", "CXCL12"],
        "immune_infiltration_thresholds": {
            "CD8_T_fraction": "0.05-0.15 (stromal)",
            "cytotoxic_score": "0.2-0.4",
            "IFNg_signature": "0.2-0.4",
            "TIL_density": "100-500_per_mm2 (peritumoral)",
        },
        "treatment_implications": [
            "Target TGF-beta pathway to promote immune penetration",
            "Anti-angiogenic therapy to normalize vasculature",
            "CAF-targeted therapy to remodel stroma",
            "Combine checkpoint inhibitor with anti-VEGF",
            "Bispecific T-cell engagers to redirect T cells",
        ],
        "clinical_trials": [
            "IMbrave150 (atezo + bev)", "MORPHEUS trials",
            "bintrafusp alfa trials",
        ],
    },
    "immunosuppressive": {
        "description": "Tumor with abundant immune infiltration dominated by "
                       "immunosuppressive cells (Tregs, MDSCs, M2 macrophages) "
                       "creating an anti-inflammatory milieu.",
        "signature_genes": ["FOXP3", "IL2RA", "TGFB1", "IL10", "ARG1",
                           "S100A9", "S100A8", "CD163", "MRC1", "IDO1",
                           "VEGFA", "CCL22", "CCL28"],
        "immune_infiltration_thresholds": {
            "Treg_fraction": ">0.10",
            "M2_macrophage_fraction": ">0.20",
            "MDSC_fraction": ">0.10",
            "CD8_to_Treg_ratio": "<2.0",
        },
        "treatment_implications": [
            "Deplete Tregs (anti-CCR4, low-dose cyclophosphamide)",
            "Reprogram M2 macrophages (anti-CSF1R, CD47-SIRPa axis)",
            "Block IDO/arginase to relieve metabolic suppression",
            "Combine checkpoint inhibitor with Treg depletion",
            "Target MDSC recruitment (anti-CXCR2, anti-CCL2)",
        ],
        "clinical_trials": [
            "ECHO-301 (epacadostat + pembro)", "TIGIT trials",
            "Magrolimab trials", "CC-90002 trials",
        ],
    },
}


# ===================================================================
# DRUG SENSITIVITY DATABASE (22 drugs)
# ===================================================================

DRUG_SENSITIVITY_DATABASE: Dict[str, Dict[str, Any]] = {
    "pembrolizumab": {
        "target": "PD-1",
        "mechanism": "Anti-PD-1 immune checkpoint inhibitor blocking the "
                     "PD-1/PD-L1 axis to restore T-cell anti-tumor activity.",
        "sensitive_cell_types": ["CD8_T", "CD4_T", "NK"],
        "resistant_cell_types": ["Treg", "M2_macrophage", "MDSC"],
        "ic50_range": "N/A (biologics)",
        "key_trials": ["KEYNOTE-024", "KEYNOTE-042", "KEYNOTE-189",
                       "KEYNOTE-590", "KEYNOTE-522"],
    },
    "nivolumab": {
        "target": "PD-1",
        "mechanism": "Anti-PD-1 monoclonal antibody restoring anti-tumor "
                     "T-cell function.",
        "sensitive_cell_types": ["CD8_T", "CD4_T", "NK"],
        "resistant_cell_types": ["Treg", "M2_macrophage"],
        "ic50_range": "N/A (biologics)",
        "key_trials": ["CheckMate 067", "CheckMate 227", "CheckMate 274",
                       "CheckMate 649"],
    },
    "ipilimumab": {
        "target": "CTLA-4",
        "mechanism": "Anti-CTLA-4 antibody enhancing T-cell activation by "
                     "blocking the inhibitory CTLA-4 receptor.",
        "sensitive_cell_types": ["CD4_T", "CD8_T"],
        "resistant_cell_types": ["Treg"],
        "ic50_range": "N/A (biologics)",
        "key_trials": ["CheckMate 067", "CA184-024", "CheckMate 214"],
    },
    "atezolizumab": {
        "target": "PD-L1",
        "mechanism": "Anti-PD-L1 antibody preventing PD-L1 engagement with "
                     "PD-1 and B7.1 on T cells.",
        "sensitive_cell_types": ["CD8_T", "CD4_T"],
        "resistant_cell_types": ["Treg", "M2_macrophage"],
        "ic50_range": "N/A (biologics)",
        "key_trials": ["IMvigor210", "IMpower150", "IMbrave150"],
    },
    "rituximab": {
        "target": "CD20",
        "mechanism": "Anti-CD20 monoclonal antibody depleting B cells via "
                     "ADCC, CDC, and direct apoptosis.",
        "sensitive_cell_types": ["B_cell", "germinal_center_B"],
        "resistant_cell_types": ["Plasma", "HSC"],
        "ic50_range": "N/A (biologics)",
        "key_trials": ["GELA-LNH98.5", "PRIMA", "GALLIUM"],
    },
    "venetoclax": {
        "target": "BCL-2",
        "mechanism": "BH3 mimetic inhibiting BCL-2 anti-apoptotic protein "
                     "to trigger mitochondrial apoptosis in cancer cells.",
        "sensitive_cell_types": ["B_cell", "CLL_cell", "AML_blast"],
        "resistant_cell_types": ["T_cell", "NK"],
        "ic50_range": "0.001-0.01 uM (CLL)",
        "key_trials": ["MURANO", "CLL14", "VIALE-A"],
    },
    "ibrutinib": {
        "target": "BTK",
        "mechanism": "Bruton's tyrosine kinase inhibitor blocking BCR "
                     "signaling in B-cell malignancies.",
        "sensitive_cell_types": ["B_cell", "CLL_cell", "MCL_cell"],
        "resistant_cell_types": ["T_cell", "myeloid"],
        "ic50_range": "0.5-5 nM",
        "key_trials": ["RESONATE", "RESONATE-2", "iLLUMINATE"],
    },
    "trastuzumab": {
        "target": "HER2",
        "mechanism": "Anti-HER2 monoclonal antibody blocking HER2 signaling "
                     "and inducing ADCC in HER2-overexpressing tumor cells.",
        "sensitive_cell_types": ["HER2_positive_epithelial"],
        "resistant_cell_types": ["HER2_negative_epithelial", "stromal"],
        "ic50_range": "N/A (biologics)",
        "key_trials": ["HERA", "CLEOPATRA", "DESTINY-Breast03"],
    },
    "erlotinib": {
        "target": "EGFR",
        "mechanism": "EGFR tyrosine kinase inhibitor blocking proliferative "
                     "signaling in EGFR-mutant tumors.",
        "sensitive_cell_types": ["EGFR_mutant_epithelial"],
        "resistant_cell_types": ["T790M_mutant", "MET_amplified"],
        "ic50_range": "2-20 nM (EGFR del19/L858R)",
        "key_trials": ["EURTAC", "OPTIMAL", "ENSURE"],
    },
    "osimertinib": {
        "target": "EGFR (T790M)",
        "mechanism": "Third-generation EGFR TKI targeting T790M resistance "
                     "mutation and sensitizing EGFR mutations.",
        "sensitive_cell_types": ["EGFR_mutant_epithelial", "T790M_mutant"],
        "resistant_cell_types": ["C797S_mutant", "MET_amplified"],
        "ic50_range": "1-15 nM",
        "key_trials": ["FLAURA", "AURA3", "ADAURA"],
    },
    "imatinib": {
        "target": "BCR-ABL / KIT / PDGFR",
        "mechanism": "Multi-kinase inhibitor targeting BCR-ABL fusion protein, "
                     "KIT, and PDGFR in CML and GIST.",
        "sensitive_cell_types": ["CML_blast", "GIST_cell"],
        "resistant_cell_types": ["T315I_mutant"],
        "ic50_range": "0.1-0.5 uM (BCR-ABL)",
        "key_trials": ["IRIS", "ENESTnd"],
    },
    "bortezomib": {
        "target": "Proteasome (26S)",
        "mechanism": "Proteasome inhibitor causing accumulation of misfolded "
                     "proteins and ER stress-mediated apoptosis in myeloma cells.",
        "sensitive_cell_types": ["Plasma", "myeloma_cell"],
        "resistant_cell_types": ["T_cell", "NK"],
        "ic50_range": "1-10 nM",
        "key_trials": ["VISTA", "ALCYONE", "CASSIOPEIA"],
    },
    "lenalidomide": {
        "target": "Cereblon (CRBN)",
        "mechanism": "IMiD promoting degradation of Ikaros/Aiolos via cereblon "
                     "E3 ligase, enhancing NK/T-cell function.",
        "sensitive_cell_types": ["myeloma_cell", "B_cell", "del5q_MDS"],
        "resistant_cell_types": ["CRBN_mutant"],
        "ic50_range": "0.1-1 uM",
        "key_trials": ["REVLIMID MM-009", "POLLUX", "MAIA"],
    },
    "azacitidine": {
        "target": "DNMT (DNA methyltransferase)",
        "mechanism": "Hypomethylating agent inducing re-expression of silenced "
                     "tumor suppressor genes and differentiation.",
        "sensitive_cell_types": ["AML_blast", "MDS_blast", "HSC"],
        "resistant_cell_types": ["mature_myeloid"],
        "ic50_range": "0.3-3 uM",
        "key_trials": ["AZA-001", "VIALE-A", "QUAZAR AML-001"],
    },
    "cytarabine": {
        "target": "DNA polymerase",
        "mechanism": "Nucleoside analog inhibiting DNA synthesis during S-phase "
                     "of cell cycle in rapidly dividing cells.",
        "sensitive_cell_types": ["AML_blast", "ALL_blast", "lymphoma_cell"],
        "resistant_cell_types": ["quiescent_HSC"],
        "ic50_range": "0.01-0.1 uM",
        "key_trials": ["MRC AML trials", "ECOG E1900", "HOVON/SAKK"],
    },
    "dexamethasone": {
        "target": "Glucocorticoid receptor",
        "mechanism": "Corticosteroid inducing apoptosis in lymphoid cells and "
                     "suppressing inflammatory cytokine production.",
        "sensitive_cell_types": ["B_cell", "T_cell", "myeloma_cell",
                                "ALL_blast"],
        "resistant_cell_types": ["myeloid", "epithelial"],
        "ic50_range": "1-100 nM (lymphoid)",
        "key_trials": ["ECOG E4A02", "MAIA", "CASSIOPEIA"],
    },
    "daratumumab": {
        "target": "CD38",
        "mechanism": "Anti-CD38 monoclonal antibody inducing tumor cell killing "
                     "via ADCC, CDC, and macrophage-mediated phagocytosis.",
        "sensitive_cell_types": ["myeloma_cell", "Plasma"],
        "resistant_cell_types": ["CD38_low_cells"],
        "ic50_range": "N/A (biologics)",
        "key_trials": ["MAIA", "CASSIOPEIA", "GRIFFIN", "APOLLO"],
    },
    "sotorasib": {
        "target": "KRAS G12C",
        "mechanism": "Covalent inhibitor locking KRAS G12C in the inactive "
                     "GDP-bound state, blocking oncogenic RAS signaling.",
        "sensitive_cell_types": ["KRAS_G12C_epithelial"],
        "resistant_cell_types": ["KRAS_wildtype", "KRAS_other_mutant"],
        "ic50_range": "10-100 nM",
        "key_trials": ["CodeBreaK 100", "CodeBreaK 200"],
    },
    "enfortumab_vedotin": {
        "target": "Nectin-4 (ADC)",
        "mechanism": "Antibody-drug conjugate targeting Nectin-4 delivering "
                     "monomethyl auristatin E (MMAE) to tumor cells.",
        "sensitive_cell_types": ["urothelial_epithelial", "Nectin4_positive"],
        "resistant_cell_types": ["Nectin4_negative"],
        "ic50_range": "0.1-1 nM (cell lines)",
        "key_trials": ["EV-301", "EV-302", "EV-103"],
    },
    "sacituzumab_govitecan": {
        "target": "Trop-2 (ADC)",
        "mechanism": "Antibody-drug conjugate targeting Trop-2 with SN-38 "
                     "payload for internalization and DNA damage.",
        "sensitive_cell_types": ["TNBC_epithelial", "urothelial_epithelial"],
        "resistant_cell_types": ["Trop2_low"],
        "ic50_range": "0.1-1 nM (cell lines)",
        "key_trials": ["ASCENT", "TROPiCS-02", "TROPHY-U-01"],
    },
    "bispecific_T_cell_engager": {
        "target": "CD3 x tumor antigen",
        "mechanism": "Bispecific antibody redirecting T cells to tumor cells "
                     "by simultaneously binding CD3 on T cells and a tumor "
                     "antigen.",
        "sensitive_cell_types": ["CD8_T", "CD4_T"],
        "resistant_cell_types": ["antigen_loss_tumor"],
        "ic50_range": "0.001-0.1 nM",
        "key_trials": ["ELREXFIO (elranatamab)", "TECVAYLI (teclistamab)",
                       "Blincyto (blinatumomab)"],
    },
    "car_t_therapy": {
        "target": "CD19 / BCMA / other",
        "mechanism": "Autologous T cells engineered with chimeric antigen "
                     "receptor to recognize and kill tumor cells expressing "
                     "the target antigen.",
        "sensitive_cell_types": ["B_cell", "myeloma_cell", "ALL_blast"],
        "resistant_cell_types": ["antigen_negative_escape",
                                "T_cell_exhaustion"],
        "ic50_range": "N/A (cell therapy)",
        "key_trials": ["ZUMA-1 (axi-cel)", "ELIANA (tisa-cel)",
                       "KarMMa (ide-cel)", "CARTITUDE-1 (cilta-cel)"],
    },
    "adagrasib": {
        "target": "KRAS G12C",
        "mechanism": "Second-generation covalent KRAS G12C inhibitor with "
                     "longer half-life and CNS penetration, locking KRAS in "
                     "the inactive GDP-bound conformation.",
        "sensitive_cell_types": ["KRAS_G12C_epithelial"],
        "resistant_cell_types": ["KRAS_wildtype", "KRAS_other_mutant",
                                "KRAS_G12C_acquired_resistance"],
        "ic50_range": "5-50 nM",
        "key_trials": ["KRYSTAL-1", "KRYSTAL-7", "KRYSTAL-12"],
    },
    "trastuzumab_deruxtecan": {
        "target": "HER2 (ADC)",
        "mechanism": "HER2-directed antibody-drug conjugate with topoisomerase "
                     "I inhibitor payload (DXd). Effective in HER2-low "
                     "expressing cells via bystander effect.",
        "sensitive_cell_types": ["HER2_positive_epithelial",
                                "HER2_low_epithelial"],
        "resistant_cell_types": ["HER2_negative_epithelial"],
        "ic50_range": "0.01-1 nM (cell lines)",
        "key_trials": ["DESTINY-Breast03", "DESTINY-Breast04",
                       "DESTINY-Lung01", "DESTINY-Gastric01"],
    },
    "belantamab_mafodotin": {
        "target": "BCMA (ADC)",
        "mechanism": "BCMA-directed antibody-drug conjugate delivering MMAF "
                     "to BCMA-expressing multiple myeloma plasma cells, "
                     "inducing mitotic arrest and apoptosis.",
        "sensitive_cell_types": ["myeloma_cell", "Plasma"],
        "resistant_cell_types": ["BCMA_low_cells"],
        "ic50_range": "0.01-0.1 nM (cell lines)",
        "key_trials": ["DREAMM-1", "DREAMM-2", "DREAMM-7", "DREAMM-8"],
    },
    "bispecific_teclistamab": {
        "target": "BCMA x CD3",
        "mechanism": "Bispecific antibody redirecting T cells to BCMA-positive "
                     "myeloma cells by simultaneously engaging CD3 on T cells "
                     "and BCMA on tumor cells.",
        "sensitive_cell_types": ["CD8_T", "CD4_T"],
        "resistant_cell_types": ["BCMA_negative_escape", "T_cell_exhaustion"],
        "ic50_range": "0.001-0.01 nM",
        "key_trials": ["MajesTEC-1", "MajesTEC-2"],
    },
    "tarlatamab": {
        "target": "DLL3 x CD3",
        "mechanism": "DLL3-targeted bispecific T-cell engager redirecting "
                     "T cells to DLL3-expressing neuroendocrine tumor cells "
                     "in small cell lung cancer.",
        "sensitive_cell_types": ["SCLC_neuroendocrine", "CD8_T"],
        "resistant_cell_types": ["DLL3_negative", "T_cell_exhaustion"],
        "ic50_range": "0.001-0.1 nM",
        "key_trials": ["DeLLphi-301", "DeLLphi-304"],
    },
    "capivasertib": {
        "target": "AKT1/2/3",
        "mechanism": "Pan-AKT inhibitor blocking PI3K/AKT/mTOR signaling in "
                     "cells with PIK3CA, AKT1, or PTEN alterations.",
        "sensitive_cell_types": ["PIK3CA_mutant_epithelial",
                                "AKT1_mutant_epithelial",
                                "PTEN_loss_epithelial"],
        "resistant_cell_types": ["KRAS_mutant", "PI3K_pathway_wildtype"],
        "ic50_range": "3-100 nM",
        "key_trials": ["CAPItello-291", "CAPItello-290"],
    },
    "elacestrant": {
        "target": "Estrogen receptor (oral SERD)",
        "mechanism": "Oral selective estrogen receptor degrader targeting "
                     "ESR1-mutant breast cancer cells resistant to standard "
                     "endocrine therapy.",
        "sensitive_cell_types": ["ER_positive_epithelial",
                                "ESR1_mutant_epithelial"],
        "resistant_cell_types": ["ER_negative_epithelial",
                                "ESR1_ligand_binding_domain_loss"],
        "ic50_range": "1-10 nM",
        "key_trials": ["EMERALD"],
    },
    "inavolisib": {
        "target": "PI3K-alpha",
        "mechanism": "Highly selective PI3K-alpha inhibitor with mutant-selective "
                     "degradation activity against PIK3CA-mutant tumor cells.",
        "sensitive_cell_types": ["PIK3CA_mutant_epithelial"],
        "resistant_cell_types": ["PIK3CA_wildtype", "PTEN_loss"],
        "ic50_range": "1-10 nM (PIK3CA-mutant)",
        "key_trials": ["INAVO120"],
    },
}


# ===================================================================
# SPATIAL TRANSCRIPTOMICS PLATFORMS
# ===================================================================

SPATIAL_PLATFORMS: Dict[str, Dict[str, Any]] = {
    "Visium": {
        "manufacturer": "10x Genomics",
        "resolution": "55 um (spot diameter, ~1-10 cells per spot)",
        "gene_count": "Whole transcriptome (~20,000 genes)",
        "coverage": "6.5 mm x 6.5 mm capture area",
        "throughput": "~5,000 spots per section",
        "analysis_methods": [
            "Spot deconvolution (cell2location, RCTD, SPOTlight)",
            "Spatially variable gene detection (SpatialDE, SPARK)",
            "Spatial clustering (BayesSpace, SpaGCN)",
            "Cell-cell communication (stLearn, COMMOT)",
            "Trajectory inference in spatial context",
        ],
        "strengths": ["Unbiased whole-transcriptome", "Established workflow",
                      "FFPE compatible (Visium CytAssist)"],
        "limitations": ["Multi-cell resolution", "Limited to tissue sections"],
    },
    "MERFISH": {
        "manufacturer": "Vizgen",
        "resolution": "100 nm (subcellular, single-molecule)",
        "gene_count": "500-1,000 genes (custom panel)",
        "coverage": "Up to 1 cm x 1 cm (MERSCOPE)",
        "throughput": "Millions of transcripts per experiment",
        "analysis_methods": [
            "Cell segmentation (Cellpose, StarDist)",
            "Spatial gene expression clustering",
            "Cell-cell interaction analysis (Squidpy)",
            "Subcellular transcript localization",
            "Spatial niche identification",
        ],
        "strengths": ["Single-molecule resolution", "High sensitivity",
                      "Subcellular localization", "Large FOV"],
        "limitations": ["Pre-selected gene panel", "Long imaging time",
                        "Computationally intensive"],
    },
    "Xenium": {
        "manufacturer": "10x Genomics",
        "resolution": "Subcellular (single-molecule FISH)",
        "gene_count": "100-5,000 genes (panel-based, customizable)",
        "coverage": "Up to 24 mm x 36 mm (multi-area)",
        "throughput": "Hundreds of thousands of cells per run",
        "analysis_methods": [
            "Cell segmentation (built-in or Cellpose)",
            "Spatial clustering and niche analysis",
            "Integration with Visium data",
            "Cell type deconvolution",
            "Spatial gene expression visualization",
        ],
        "strengths": ["Subcellular resolution", "FFPE compatible",
                      "Expanding panel sizes", "Integrated workflow"],
        "limitations": ["Panel-based (not whole transcriptome)",
                        "Higher cost per gene than Visium"],
    },
    "CODEX": {
        "manufacturer": "Akoya Biosciences",
        "resolution": "Single-cell (protein-level, ~250 nm optical)",
        "gene_count": "40-100 protein markers (antibody panels)",
        "coverage": "Full tissue section",
        "throughput": "Millions of cells per tissue section",
        "analysis_methods": [
            "Cell segmentation (Mesmer, DeepCell)",
            "Phenotyping and clustering (Phenograph, FlowSOM)",
            "Spatial neighborhood analysis",
            "Cell-cell proximity and interaction scoring",
            "Tissue microenvironment classification",
        ],
        "strengths": ["Protein-level detection", "High multiplexing",
                      "Compatible with FFPE", "Phenotypic resolution"],
        "limitations": ["Antibody availability limits targets",
                        "No transcript-level information",
                        "Requires antibody validation"],
    },
}


# ===================================================================
# MARKER GENE DATABASE (55 canonical markers)
# ===================================================================

MARKER_GENE_DATABASE: Dict[str, Dict[str, Any]] = {
    "CD3D": {"cell_types": ["T_cell", "CD4_T", "CD8_T", "Treg", "NKT"],
             "function": "T-cell receptor complex component"},
    "CD3E": {"cell_types": ["T_cell", "CD4_T", "CD8_T", "Treg"],
             "function": "T-cell receptor signaling subunit"},
    "CD8A": {"cell_types": ["CD8_T"],
             "function": "MHC class I coreceptor"},
    "CD8B": {"cell_types": ["CD8_T"],
             "function": "MHC class I coreceptor beta chain"},
    "CD4": {"cell_types": ["CD4_T", "Treg", "Monocyte"],
            "function": "MHC class II coreceptor and HIV receptor"},
    "FOXP3": {"cell_types": ["Treg"],
              "function": "Master transcription factor for regulatory T cells"},
    "IL2RA": {"cell_types": ["Treg", "activated_T"],
              "function": "IL-2 receptor alpha chain (CD25)"},
    "CD19": {"cell_types": ["B_cell"],
             "function": "B-cell surface marker and BCR coreceptor"},
    "MS4A1": {"cell_types": ["B_cell"],
              "function": "CD20, target of rituximab"},
    "CD79A": {"cell_types": ["B_cell"],
              "function": "B-cell receptor signaling component"},
    "SDC1": {"cell_types": ["Plasma"],
             "function": "CD138, plasma cell surface proteoglycan"},
    "MZB1": {"cell_types": ["Plasma", "B_cell"],
             "function": "ER-resident chaperone in antibody-secreting cells"},
    "KLRD1": {"cell_types": ["NK", "CD8_T"],
              "function": "CD94, NK cell receptor component"},
    "NKG7": {"cell_types": ["NK", "CD8_T", "gamma_delta_T"],
             "function": "Cytotoxic granule membrane protein"},
    "NCAM1": {"cell_types": ["NK"],
              "function": "CD56, neural cell adhesion molecule on NK cells"},
    "CD14": {"cell_types": ["Monocyte", "Macrophage"],
             "function": "LPS coreceptor on myeloid cells"},
    "LYZ": {"cell_types": ["Monocyte", "Macrophage", "Neutrophil"],
            "function": "Lysozyme, antimicrobial enzyme"},
    "CD68": {"cell_types": ["Macrophage"],
             "function": "Macrophage-associated glycoprotein"},
    "MARCO": {"cell_types": ["Macrophage"],
              "function": "Macrophage scavenger receptor"},
    "CSF1R": {"cell_types": ["Macrophage", "Monocyte"],
              "function": "M-CSF receptor driving macrophage differentiation"},
    "CLEC9A": {"cell_types": ["cDC1"],
               "function": "C-type lectin on cross-presenting DCs"},
    "CD1C": {"cell_types": ["cDC2"],
             "function": "Lipid antigen presentation molecule"},
    "CSF3R": {"cell_types": ["Neutrophil"],
              "function": "G-CSF receptor controlling neutrophil production"},
    "S100A8": {"cell_types": ["Neutrophil", "Monocyte", "MDSC"],
               "function": "Calprotectin subunit, alarmin"},
    "S100A9": {"cell_types": ["Neutrophil", "Monocyte", "MDSC"],
               "function": "Calprotectin subunit, alarmin"},
    "COL1A1": {"cell_types": ["Fibroblast"],
               "function": "Type I collagen alpha-1 chain"},
    "DCN": {"cell_types": ["Fibroblast"],
            "function": "Decorin, collagen-binding proteoglycan"},
    "FAP": {"cell_types": ["cancer_associated_fibroblast"],
            "function": "Fibroblast activation protein"},
    "PECAM1": {"cell_types": ["Endothelial"],
               "function": "CD31, platelet-endothelial adhesion molecule"},
    "VWF": {"cell_types": ["Endothelial"],
            "function": "Von Willebrand factor, hemostasis mediator"},
    "CDH5": {"cell_types": ["Endothelial"],
             "function": "VE-cadherin, endothelial adherens junction"},
    "EPCAM": {"cell_types": ["Epithelial"],
              "function": "Epithelial cell adhesion molecule"},
    "KRT18": {"cell_types": ["Epithelial", "Hepatocyte"],
              "function": "Cytokeratin 18, epithelial intermediate filament"},
    "KRT19": {"cell_types": ["Epithelial"],
              "function": "Cytokeratin 19, ductal/luminal epithelial marker"},
    "ALB": {"cell_types": ["Hepatocyte"],
            "function": "Albumin, major serum protein"},
    "APOB": {"cell_types": ["Hepatocyte"],
             "function": "Apolipoprotein B, lipoprotein component"},
    "MAP2": {"cell_types": ["Neuron"],
             "function": "Microtubule-associated protein 2, neuronal marker"},
    "RBFOX3": {"cell_types": ["Neuron"],
               "function": "NeuN, post-mitotic neuron marker"},
    "GFAP": {"cell_types": ["Astrocyte"],
             "function": "Glial fibrillary acidic protein"},
    "AQP4": {"cell_types": ["Astrocyte"],
             "function": "Aquaporin-4, water channel in astrocytic endfeet"},
    "MBP": {"cell_types": ["Oligodendrocyte"],
            "function": "Myelin basic protein"},
    "PLP1": {"cell_types": ["Oligodendrocyte"],
             "function": "Proteolipid protein 1, major myelin component"},
    "TNNT2": {"cell_types": ["Cardiomyocyte"],
              "function": "Cardiac troponin T"},
    "MYH6": {"cell_types": ["Cardiomyocyte"],
             "function": "Myosin heavy chain 6 (alpha)"},
    "ACTA2": {"cell_types": ["Smooth_muscle", "myofibroblast", "Pericyte"],
              "function": "Alpha-smooth muscle actin"},
    "MYH11": {"cell_types": ["Smooth_muscle"],
              "function": "Smooth muscle myosin heavy chain"},
    "ADIPOQ": {"cell_types": ["Adipocyte"],
               "function": "Adiponectin, adipokine"},
    "FABP4": {"cell_types": ["Adipocyte"],
              "function": "Fatty acid binding protein 4"},
    "CD34": {"cell_types": ["HSC", "Endothelial"],
             "function": "Hematopoietic stem cell and endothelial marker"},
    "KIT": {"cell_types": ["HSC", "Mast_cell"],
            "function": "CD117, stem cell factor receptor"},
    "PDGFRB": {"cell_types": ["Pericyte", "Fibroblast"],
               "function": "PDGF receptor beta, mural cell marker"},
    "CX3CR1": {"cell_types": ["Microglia", "Monocyte"],
               "function": "Fractalkine receptor"},
    "P2RY12": {"cell_types": ["Microglia"],
               "function": "Purinergic receptor, homeostatic microglia marker"},
    "TMEM119": {"cell_types": ["Microglia"],
                "function": "Microglia-specific transmembrane protein"},
    "HBA1": {"cell_types": ["Erythrocyte_precursor"],
             "function": "Hemoglobin alpha 1"},
    "TRGV9": {"cell_types": ["Gamma_delta_T"],
              "function": "T-cell receptor gamma variable 9, Vgamma9Vdelta2 "
                          "gamma-delta T-cell marker"},
    "SLC4A10": {"cell_types": ["MAIT_cell"],
                "function": "Sodium bicarbonate transporter, MAIT cell marker"},
    "CLEC4C": {"cell_types": ["Plasmacytoid_DC", "pDC"],
               "function": "C-type lectin receptor on plasmacytoid DCs (BDCA-2)"},
    "XCR1": {"cell_types": ["cDC1"],
             "function": "Chemokine receptor on cross-presenting cDC1 cells"},
    "RORC": {"cell_types": ["Innate_lymphoid_cell_3", "Th17"],
             "function": "RAR-related orphan receptor C, master TF for ILC3/Th17"},
    "GATA3_ILC": {"cell_types": ["Innate_lymphoid_cell_2", "Th2"],
                  "function": "GATA binding protein 3 in ILC2 context (also Th2)"},
    "ITGA2B": {"cell_types": ["Megakaryocyte"],
               "function": "Integrin alpha-2b (CD41), platelet glycoprotein IIb"},
    "GYPA": {"cell_types": ["Erythroid_progenitor", "Erythrocyte_precursor"],
             "function": "Glycophorin A (CD235a), erythroid lineage marker"},
    "RGS5": {"cell_types": ["Pericyte"],
             "function": "Regulator of G-protein signaling 5, pericyte marker"},
    "MPZ": {"cell_types": ["Schwann_cell"],
            "function": "Myelin protein zero, major PNS myelin component"},
    "NPHS1": {"cell_types": ["Podocyte"],
              "function": "Nephrin, slit diaphragm component essential for "
                          "glomerular filtration"},
    "MUC2": {"cell_types": ["Goblet_cell"],
             "function": "Mucin 2, gel-forming mucin of intestinal mucus layer"},
    "POSTN": {"cell_types": ["cancer_associated_fibroblast", "Fibroblast"],
              "function": "Periostin, ECM protein promoting tumor invasion and "
                          "stromal remodeling"},
    "TCF7": {"cell_types": ["stem_memory_T", "naive_T", "progenitor_exhausted_T"],
             "function": "TCF1, transcription factor marking stem-like memory "
                         "T cells"},
    "CXCR6": {"cell_types": ["tissue_resident_memory_T", "NKT"],
              "function": "Chemokine receptor for CXCL16, tissue-resident "
                          "memory T-cell marker"},
    "MKI67": {"cell_types": ["proliferating_cell"],
              "function": "Ki-67, universal marker of cell proliferation "
                          "(absent in G0)"},
    "TOP2A": {"cell_types": ["proliferating_cell"],
              "function": "Topoisomerase II alpha, S/G2/M phase proliferation "
                          "marker"},
    "CDK1": {"cell_types": ["proliferating_cell"],
             "function": "Cyclin-dependent kinase 1, master regulator of "
                         "mitotic entry"},
    "PCNA": {"cell_types": ["proliferating_cell"],
             "function": "Proliferating cell nuclear antigen, DNA replication "
                         "processivity factor"},
    "MCM2": {"cell_types": ["proliferating_cell"],
             "function": "Minichromosome maintenance complex component 2, "
                         "replication licensing factor"},
}


# ===================================================================
# IMMUNE SIGNATURES
# ===================================================================

IMMUNE_SIGNATURES: Dict[str, Dict[str, Any]] = {
    "cytotoxic": {
        "genes": ["GZMA", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7",
                  "GNLY", "FASLG"],
        "description": "Cytotoxic effector signature associated with active "
                       "anti-tumor killing via perforin/granzyme pathway.",
        "cell_types": ["CD8_T", "NK", "gamma_delta_T"],
        "clinical_relevance": "High cytotoxic score correlates with checkpoint "
                              "inhibitor response and improved OS in multiple "
                              "tumor types.",
    },
    "exhaustion": {
        "genes": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX",
                  "ENTPD1", "BATF"],
        "description": "T-cell exhaustion signature with co-expression of "
                       "multiple inhibitory receptors and transcription factor "
                       "TOX.",
        "cell_types": ["CD8_T", "CD4_T"],
        "clinical_relevance": "Exhausted T cells retain partial function and "
                              "are the primary target of checkpoint inhibitor "
                              "therapy. PD-1+TCF1+ progenitor exhausted cells "
                              "predict ICB response.",
    },
    "regulatory": {
        "genes": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2",
                  "IL10", "TGFB1"],
        "description": "Regulatory T-cell signature suppressing anti-tumor "
                       "immunity through contact-dependent and cytokine-mediated "
                       "mechanisms.",
        "cell_types": ["Treg"],
        "clinical_relevance": "High Treg infiltration in tumor associates with "
                              "poor prognosis in most solid tumors. High "
                              "CD8/Treg ratio is a favorable prognostic marker.",
    },
    "myeloid_suppression": {
        "genes": ["S100A9", "S100A8", "ARG1", "ARG2", "IDO1", "NOS2",
                  "CD163", "MRC1"],
        "description": "Myeloid-derived suppressive signature from MDSCs and "
                       "M2-polarized macrophages creating an immunosuppressive "
                       "tumor microenvironment.",
        "cell_types": ["MDSC", "M2_macrophage",
                       "tumor_associated_macrophage"],
        "clinical_relevance": "Myeloid suppression correlates with ICB "
                              "resistance. Anti-CSF1R and ARG inhibitors are "
                              "being tested to reprogram suppressive myeloid "
                              "compartment.",
    },
    "fibroblast_activation": {
        "genes": ["FAP", "ACTA2", "COL1A1", "COL1A2", "POSTN", "FN1",
                  "TGFB1", "TGFB2"],
        "description": "Cancer-associated fibroblast activation signature "
                       "indicating stromal remodeling and immune exclusion.",
        "cell_types": ["cancer_associated_fibroblast", "myofibroblast"],
        "clinical_relevance": "High CAF signature predicts immune exclusion "
                              "phenotype and poor response to immunotherapy. "
                              "TGF-beta blockade may overcome stromal barrier.",
    },
    "memory_T": {
        "genes": ["TCF7", "IL7R", "CCR7", "SELL", "LEF1"],
        "description": "Stem-like memory T-cell signature associated with "
                       "durable anti-tumor responses and long-term immunity. "
                       "TCF1+ progenitor cells sustain effector T-cell pools.",
        "cell_types": ["memory_CD8", "memory_CD4", "stem_memory_T"],
        "clinical_relevance": "High TCF7+/PD-1+ progenitor exhausted cells "
                              "predict durable response to checkpoint inhibitor "
                              "therapy. Associated with improved OS across "
                              "multiple tumor types.",
    },
    "tissue_resident_memory": {
        "genes": ["ITGAE", "ZNF683", "CXCR6", "CD69", "ITGA1"],
        "description": "Tissue-resident memory T-cell signature marking "
                       "non-circulating T cells that provide local tissue "
                       "surveillance without re-entering circulation.",
        "cell_types": ["tissue_resident_memory_CD8", "tissue_resident_memory_CD4"],
        "clinical_relevance": "CD103+ tissue-resident memory T cells in tumor "
                              "correlate with improved prognosis and response to "
                              "immunotherapy in NSCLC, melanoma, and breast cancer.",
    },
    "M1_macrophage": {
        "genes": ["NOS2", "TNF", "IL1B", "IL6", "CD80"],
        "description": "Pro-inflammatory M1 macrophage polarization signature "
                       "associated with anti-tumor immunity, pathogen clearance, "
                       "and Th1 immune responses.",
        "cell_types": ["M1_macrophage", "activated_macrophage"],
        "clinical_relevance": "High M1/M2 ratio in tumor associates with better "
                              "prognosis and improved checkpoint inhibitor response. "
                              "M1 repolarization is a therapeutic strategy.",
    },
    "M2_macrophage": {
        "genes": ["MRC1", "CD163", "TGFB1", "IL10", "CCL18"],
        "description": "Anti-inflammatory M2 macrophage polarization signature "
                       "associated with tissue remodeling, immunosuppression, "
                       "and pro-tumorigenic functions.",
        "cell_types": ["M2_macrophage", "tumor_associated_macrophage"],
        "clinical_relevance": "M2-dominant tumor macrophages predict poor prognosis "
                              "and ICB resistance. Anti-CSF1R, CD47/SIRPa blockade, "
                              "and PI3Kg inhibition aim to reprogram M2 to M1.",
    },
    "cancer_associated_fibroblast": {
        "genes": ["FAP", "ACTA2", "PDGFRA", "COL1A1", "POSTN"],
        "description": "Cancer-associated fibroblast signature marking stromal "
                       "cells that create a physical and immunological barrier "
                       "preventing immune cell infiltration.",
        "cell_types": ["cancer_associated_fibroblast", "myofibroblast"],
        "clinical_relevance": "High CAF abundance creates immune-excluded phenotype "
                              "resistant to immunotherapy. FAP-targeted therapies "
                              "and TGF-beta blockade are under investigation.",
    },
}


# ===================================================================
# FOUNDATION MODELS
# ===================================================================

FOUNDATION_MODELS: Dict[str, Dict[str, Any]] = {
    "scGPT": {
        "full_name": "Single-Cell Generative Pre-trained Transformer",
        "parameters": "~50M",
        "training_data": "33M human single-cell transcriptomes from CellXGene",
        "architecture": "Transformer with gene-token embedding and attention "
                        "masking",
        "capabilities": [
            "Cell type annotation",
            "Gene perturbation prediction",
            "Multi-batch integration",
            "Multi-omic integration",
            "Gene regulatory network inference",
            "Cell embedding generation",
        ],
        "input_format": "Gene expression vectors (top variable genes)",
        "output_format": "Cell embeddings, gene embeddings, predicted "
                         "expression",
        "reference": "Cui et al., Nature Methods 2024",
        "gpu_requirements": "1x A100 40GB for inference, 4x A100 80GB for "
                            "fine-tuning",
    },
    "Geneformer": {
        "full_name": "Geneformer",
        "parameters": "~10M",
        "training_data": "~30M human single-cell transcriptomes from "
                         "Genecorpus-30M",
        "architecture": "BERT-style transformer with rank-value encoding of "
                        "genes",
        "capabilities": [
            "Cell type classification",
            "Gene dosage sensitivity prediction",
            "Chromatin dynamics prediction",
            "In silico perturbation",
            "Disease state classification",
            "Network biology insights",
        ],
        "input_format": "Rank-ordered gene tokens per cell",
        "output_format": "Cell embeddings, gene embeddings, classification "
                         "logits",
        "reference": "Theodoris et al., Nature 2023",
        "gpu_requirements": "1x V100 32GB for inference, 2x A100 40GB for "
                            "fine-tuning",
    },
    "scFoundation": {
        "full_name": "Large-Scale Foundation Model for Single-Cell "
                     "Transcriptomics",
        "parameters": "~100M",
        "training_data": "50M+ human single-cell transcriptomes",
        "architecture": "Asymmetric transformer encoder-decoder with gene "
                        "expression value embedding",
        "capabilities": [
            "Cell type annotation across tissues",
            "Gene expression imputation",
            "Drug response prediction",
            "Perturbation response modeling",
            "Cross-tissue cell embedding",
            "Zero-shot cell type transfer",
        ],
        "input_format": "Binned gene expression values with gene tokens",
        "output_format": "Cell embeddings, imputed expression, prediction "
                         "scores",
        "reference": "Hao et al., Nature Methods 2024",
        "gpu_requirements": "1x A100 40GB for inference, 8x A100 80GB for "
                            "fine-tuning",
    },
}


# ===================================================================
# GPU BENCHMARKS
# ===================================================================

GPU_BENCHMARKS: Dict[str, Dict[str, Any]] = {
    "10K_cells": {
        "cell_count": 10_000,
        "preprocessing": {
            "normalization_sec": 0.5,
            "hvg_selection_sec": 1.2,
            "pca_sec": 0.8,
            "total_sec": 2.5,
        },
        "clustering": {
            "leiden_sec": 1.0,
            "umap_sec": 2.0,
            "total_sec": 3.0,
        },
        "differential_expression": {
            "wilcoxon_sec": 3.0,
            "total_sec": 3.0,
        },
        "foundation_model_inference": {
            "scGPT_sec": 5.0,
            "Geneformer_sec": 3.0,
            "total_sec": 8.0,
        },
        "total_pipeline_sec": 16.5,
        "gpu_memory_gb": 2.0,
        "gpu": "NVIDIA A100 40GB",
    },
    "50K_cells": {
        "cell_count": 50_000,
        "preprocessing": {
            "normalization_sec": 2.0,
            "hvg_selection_sec": 4.5,
            "pca_sec": 3.5,
            "total_sec": 10.0,
        },
        "clustering": {
            "leiden_sec": 5.0,
            "umap_sec": 12.0,
            "total_sec": 17.0,
        },
        "differential_expression": {
            "wilcoxon_sec": 15.0,
            "total_sec": 15.0,
        },
        "foundation_model_inference": {
            "scGPT_sec": 25.0,
            "Geneformer_sec": 15.0,
            "total_sec": 40.0,
        },
        "total_pipeline_sec": 82.0,
        "gpu_memory_gb": 6.0,
        "gpu": "NVIDIA A100 40GB",
    },
    "100K_cells": {
        "cell_count": 100_000,
        "preprocessing": {
            "normalization_sec": 4.0,
            "hvg_selection_sec": 9.0,
            "pca_sec": 7.0,
            "total_sec": 20.0,
        },
        "clustering": {
            "leiden_sec": 12.0,
            "umap_sec": 30.0,
            "total_sec": 42.0,
        },
        "differential_expression": {
            "wilcoxon_sec": 35.0,
            "total_sec": 35.0,
        },
        "foundation_model_inference": {
            "scGPT_sec": 55.0,
            "Geneformer_sec": 30.0,
            "total_sec": 85.0,
        },
        "total_pipeline_sec": 182.0,
        "gpu_memory_gb": 12.0,
        "gpu": "NVIDIA A100 40GB",
    },
    "500K_cells": {
        "cell_count": 500_000,
        "preprocessing": {
            "normalization_sec": 20.0,
            "hvg_selection_sec": 45.0,
            "pca_sec": 35.0,
            "total_sec": 100.0,
        },
        "clustering": {
            "leiden_sec": 60.0,
            "umap_sec": 180.0,
            "total_sec": 240.0,
        },
        "differential_expression": {
            "wilcoxon_sec": 200.0,
            "total_sec": 200.0,
        },
        "foundation_model_inference": {
            "scGPT_sec": 300.0,
            "Geneformer_sec": 160.0,
            "total_sec": 460.0,
        },
        "total_pipeline_sec": 1000.0,
        "gpu_memory_gb": 38.0,
        "gpu": "NVIDIA A100 40GB",
    },
}


# ===================================================================
# LIGAND-RECEPTOR PAIRS (25 pairs)
# ===================================================================

LIGAND_RECEPTOR_PAIRS: Dict[str, Dict[str, Any]] = {
    "CXCL12_CXCR4": {
        "ligand": "CXCL12",
        "receptor": "CXCR4",
        "function": "Stromal-hematopoietic axis regulating HSC retention in "
                    "bone marrow niche and immune cell trafficking.",
        "pathway": "Chemokine signaling",
    },
    "CCL19_CCR7": {
        "ligand": "CCL19",
        "receptor": "CCR7",
        "function": "Lymph node homing signal directing naive T cells and "
                    "mature dendritic cells to T-cell zones.",
        "pathway": "Chemokine signaling",
    },
    "PDCD1_CD274": {
        "ligand": "CD274",
        "receptor": "PDCD1",
        "function": "PD-1/PD-L1 immune checkpoint axis suppressing T-cell "
                    "effector function in tumors.",
        "pathway": "Immune checkpoint",
    },
    "CTLA4_CD80": {
        "ligand": "CD80",
        "receptor": "CTLA4",
        "function": "Co-inhibitory interaction outcompeting CD28 for B7 "
                    "ligand binding to suppress T-cell activation.",
        "pathway": "Immune checkpoint",
    },
    "CD28_CD80": {
        "ligand": "CD80",
        "receptor": "CD28",
        "function": "Co-stimulatory signal providing second activation signal "
                    "for naive T-cell priming.",
        "pathway": "Co-stimulatory",
    },
    "TIGIT_PVR": {
        "ligand": "PVR",
        "receptor": "TIGIT",
        "function": "Checkpoint interaction via CD155 suppressing NK and "
                    "T-cell cytotoxicity in the tumor microenvironment.",
        "pathway": "Immune checkpoint",
    },
    "CD40_CD40LG": {
        "ligand": "CD40LG",
        "receptor": "CD40",
        "function": "B-cell activation and germinal center formation signal "
                    "from T-helper cells to B cells and APCs.",
        "pathway": "Co-stimulatory",
    },
    "VEGFA_FLT1": {
        "ligand": "VEGFA",
        "receptor": "FLT1",
        "function": "VEGFR1-mediated angiogenesis and vascular permeability "
                    "signaling in tumor neovascularization.",
        "pathway": "Angiogenesis",
    },
    "VEGFA_KDR": {
        "ligand": "VEGFA",
        "receptor": "KDR",
        "function": "VEGFR2-mediated primary angiogenic signaling driving "
                    "endothelial cell proliferation and migration.",
        "pathway": "Angiogenesis",
    },
    "TGFB1_TGFBR2": {
        "ligand": "TGFB1",
        "receptor": "TGFBR2",
        "function": "TGF-beta immunosuppression and epithelial-mesenchymal "
                    "transition signaling in tumor stroma.",
        "pathway": "TGF-beta signaling",
    },
    "TNF_TNFRSF1A": {
        "ligand": "TNF",
        "receptor": "TNFRSF1A",
        "function": "Pro-inflammatory TNF signaling via TNFR1 mediating "
                    "apoptosis, NF-kB activation, and inflammation.",
        "pathway": "Inflammatory signaling",
    },
    "IL6_IL6R": {
        "ligand": "IL6",
        "receptor": "IL6R",
        "function": "IL-6 inflammatory and acute phase response signaling "
                    "via classical and trans-signaling pathways.",
        "pathway": "Inflammatory signaling",
    },
    "IFNG_IFNGR1": {
        "ligand": "IFNG",
        "receptor": "IFNGR1",
        "function": "IFN-gamma immune activation inducing MHC upregulation, "
                    "macrophage activation, and anti-tumor immunity.",
        "pathway": "Interferon signaling",
    },
    "IL10_IL10RA": {
        "ligand": "IL10",
        "receptor": "IL10RA",
        "function": "Anti-inflammatory immunosuppressive signaling from "
                    "Tregs and M2 macrophages dampening immune responses.",
        "pathway": "Immunosuppression",
    },
    "CSF1_CSF1R": {
        "ligand": "CSF1",
        "receptor": "CSF1R",
        "function": "Macrophage recruitment, differentiation, and survival "
                    "signaling in the tumor microenvironment.",
        "pathway": "Myeloid signaling",
    },
    "CCL2_CCR2": {
        "ligand": "CCL2",
        "receptor": "CCR2",
        "function": "Monocyte and macrophage recruitment to sites of "
                    "inflammation and tumor tissues.",
        "pathway": "Chemokine signaling",
    },
    "CXCL9_CXCR3": {
        "ligand": "CXCL9",
        "receptor": "CXCR3",
        "function": "IFN-gamma-induced T-cell chemotaxis to inflamed and "
                    "tumor tissues.",
        "pathway": "Chemokine signaling",
    },
    "CXCL10_CXCR3": {
        "ligand": "CXCL10",
        "receptor": "CXCR3",
        "function": "T-cell and NK cell chemotaxis to sites of viral "
                    "infection and tumor inflammation.",
        "pathway": "Chemokine signaling",
    },
    "FAS_FASLG": {
        "ligand": "FASLG",
        "receptor": "FAS",
        "function": "Extrinsic apoptosis pathway triggering caspase cascade "
                    "in target cells upon death receptor engagement.",
        "pathway": "Apoptosis",
    },
    "TRAIL_TNFRSF10A": {
        "ligand": "TNFSF10",
        "receptor": "TNFRSF10A",
        "function": "TRAIL/DR4 apoptosis signaling selectively killing "
                    "tumor cells while sparing normal cells.",
        "pathway": "Apoptosis",
    },
    "HGF_MET": {
        "ligand": "HGF",
        "receptor": "MET",
        "function": "Hepatocyte growth factor signaling promoting cell "
                    "growth, motility, and invasion in cancer.",
        "pathway": "Growth factor signaling",
    },
    "WNT5A_FZD5": {
        "ligand": "WNT5A",
        "receptor": "FZD5",
        "function": "Non-canonical Wnt signaling regulating cell polarity, "
                    "migration, and tissue patterning.",
        "pathway": "Wnt signaling",
    },
    "NOTCH1_DLL1": {
        "ligand": "DLL1",
        "receptor": "NOTCH1",
        "function": "Notch signaling controlling cell fate decisions, "
                    "differentiation, and stem cell maintenance.",
        "pathway": "Notch signaling",
    },
    "EGF_EGFR": {
        "ligand": "EGF",
        "receptor": "EGFR",
        "function": "Epidermal growth factor signaling driving epithelial "
                    "cell proliferation and survival.",
        "pathway": "Growth factor signaling",
    },
    "BMP2_BMPR1A": {
        "ligand": "BMP2",
        "receptor": "BMPR1A",
        "function": "BMP signaling for osteogenic differentiation, tissue "
                    "patterning, and stem cell regulation.",
        "pathway": "BMP/TGF-beta superfamily",
    },
}


# ===================================================================
# CANCER TME ATLAS (12 cancer types)
# ===================================================================

CANCER_TME_ATLAS: Dict[str, Dict[str, Any]] = {
    "NSCLC": {
        "dominant_tme_class": "variable (hot/cold/excluded)",
        "key_immune_features": [
            "High PD-L1 expression in adenocarcinoma subtype",
            "Tertiary lymphoid structures predict ICB response",
            "Smoking-associated high TMB enhances neoantigen load",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 25, "B_cells_pct": 8, "NK_pct": 5,
            "myeloid_pct": 35, "fibroblast_pct": 15,
        },
        "checkpoint_expression": ["PD-L1", "CTLA-4", "TIGIT", "LAG-3"],
        "treatment_response_pattern": "PD-L1 TPS >= 50% predicts anti-PD-1 "
                                      "monotherapy benefit; TMB-high is "
                                      "independent biomarker.",
    },
    "breast": {
        "dominant_tme_class": "subtype-dependent (TNBC hot, ER+ cold)",
        "key_immune_features": [
            "TNBC has highest TIL infiltration among subtypes",
            "ER+ tumors are generally immune-cold",
            "HER2+ shows intermediate immune infiltration",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 15, "B_cells_pct": 5, "NK_pct": 3,
            "myeloid_pct": 30, "fibroblast_pct": 25,
        },
        "checkpoint_expression": ["PD-L1 (TNBC)", "CTLA-4"],
        "treatment_response_pattern": "ICB effective in PD-L1+ TNBC with "
                                      "chemotherapy combination; limited "
                                      "benefit in ER+ disease.",
    },
    "colorectal": {
        "dominant_tme_class": "MSI-H hot / MSS cold",
        "key_immune_features": [
            "MSI-H/dMMR tumors have high neoantigen load and dense TILs",
            "MSS tumors are generally immune-cold",
            "Right-sided tumors more likely MSI-H",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 20, "B_cells_pct": 5, "NK_pct": 3,
            "myeloid_pct": 30, "fibroblast_pct": 20,
        },
        "checkpoint_expression": ["PD-1", "PD-L1 (MSI-H)", "CTLA-4"],
        "treatment_response_pattern": "MSI-H/dMMR highly responsive to "
                                      "anti-PD-1; MSS resistant to single-agent "
                                      "ICB.",
    },
    "melanoma": {
        "dominant_tme_class": "hot",
        "key_immune_features": [
            "Highest TMB among solid tumors (UV mutagenesis)",
            "Dense CD8+ T-cell infiltration",
            "Frequent PD-L1 expression and IFN-gamma signature",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 30, "B_cells_pct": 5, "NK_pct": 5,
            "myeloid_pct": 25, "fibroblast_pct": 10,
        },
        "checkpoint_expression": ["PD-1", "PD-L1", "CTLA-4", "LAG-3", "TIGIT"],
        "treatment_response_pattern": "Dual checkpoint blockade (nivo+ipi) "
                                      "achieves ~60% response rate; IFN-gamma "
                                      "signature predicts benefit.",
    },
    "PDAC": {
        "dominant_tme_class": "excluded / immunosuppressive",
        "key_immune_features": [
            "Dense desmoplastic stroma excludes T cells",
            "High CAF and M2 macrophage infiltration",
            "Low TMB and minimal neoantigen load",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 5, "B_cells_pct": 2, "NK_pct": 1,
            "myeloid_pct": 25, "fibroblast_pct": 50,
        },
        "checkpoint_expression": ["PD-L1 (low)", "CTLA-4 (low)"],
        "treatment_response_pattern": "Refractory to ICB monotherapy; requires "
                                      "stromal remodeling and combination "
                                      "strategies.",
    },
    "GBM": {
        "dominant_tme_class": "immunosuppressive",
        "key_immune_features": [
            "Blood-brain barrier limits immune infiltration",
            "High M2 macrophage/microglia dominance",
            "T-cell exhaustion and lymphopenia",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 5, "B_cells_pct": 1, "NK_pct": 1,
            "myeloid_pct": 50, "fibroblast_pct": 5,
        },
        "checkpoint_expression": ["PD-L1 (variable)", "TIM-3", "IDO1"],
        "treatment_response_pattern": "ICB has failed in unselected GBM; "
                                      "dMMR subset may benefit. CAR-T and "
                                      "oncolytic virus approaches under study.",
    },
    "HCC": {
        "dominant_tme_class": "variable (excluded to immunosuppressive)",
        "key_immune_features": [
            "Chronic inflammation background (HBV/HCV/NASH)",
            "Tolerogenic liver environment with Kupffer cells",
            "High Treg and exhausted T-cell infiltration",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 15, "B_cells_pct": 3, "NK_pct": 8,
            "myeloid_pct": 35, "fibroblast_pct": 15,
        },
        "checkpoint_expression": ["PD-1", "PD-L1", "CTLA-4", "TIM-3"],
        "treatment_response_pattern": "Atezolizumab + bevacizumab is first-line "
                                      "standard; anti-VEGF normalizes vasculature "
                                      "to enable T-cell infiltration.",
    },
    "RCC": {
        "dominant_tme_class": "hot / immunosuppressive",
        "key_immune_features": [
            "High immune infiltration but with suppressive myeloid cells",
            "VHL loss drives VEGF/HIF pathway and angiogenesis",
            "Sarcomatoid features predict ICB benefit",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 25, "B_cells_pct": 5, "NK_pct": 5,
            "myeloid_pct": 30, "fibroblast_pct": 10,
        },
        "checkpoint_expression": ["PD-1", "PD-L1", "CTLA-4", "LAG-3"],
        "treatment_response_pattern": "ICB + anti-VEGF TKI is standard "
                                      "first-line; dual ICB (nivo+ipi) for "
                                      "intermediate/poor risk.",
    },
    "ovarian": {
        "dominant_tme_class": "excluded / immunosuppressive",
        "key_immune_features": [
            "Peritoneal dissemination with ascites TME",
            "BRCA-mutant tumors have higher TIL infiltration",
            "High CAF and mesothelial cell abundance",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 10, "B_cells_pct": 5, "NK_pct": 3,
            "myeloid_pct": 30, "fibroblast_pct": 25,
        },
        "checkpoint_expression": ["PD-L1 (variable)", "CTLA-4"],
        "treatment_response_pattern": "Limited ICB benefit in unselected "
                                      "population; BRCA-mutant and high-TIL "
                                      "subsets may respond. PARP + ICB combos "
                                      "under investigation.",
    },
    "HNSCC": {
        "dominant_tme_class": "hot (HPV+) / cold (HPV-)",
        "key_immune_features": [
            "HPV-positive tumors have higher immune infiltration",
            "HPV-negative tumors have smoking-associated TMB",
            "Tertiary lymphoid structures correlate with response",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 20, "B_cells_pct": 8, "NK_pct": 5,
            "myeloid_pct": 30, "fibroblast_pct": 15,
        },
        "checkpoint_expression": ["PD-1", "PD-L1", "CTLA-4"],
        "treatment_response_pattern": "Pembrolizumab first-line for PD-L1 CPS "
                                      ">= 1; HPV+ tumors have better prognosis "
                                      "but similar ICB response rate.",
    },
    "bladder": {
        "dominant_tme_class": "variable",
        "key_immune_features": [
            "High TMB from APOBEC mutagenesis",
            "Luminal-infiltrated subtype predicts ICB benefit",
            "Nectin-4 and Trop-2 expression for ADC targeting",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 20, "B_cells_pct": 5, "NK_pct": 3,
            "myeloid_pct": 30, "fibroblast_pct": 20,
        },
        "checkpoint_expression": ["PD-1", "PD-L1", "CTLA-4"],
        "treatment_response_pattern": "ICB second-line standard; ADCs "
                                      "(enfortumab vedotin) now first-line "
                                      "with pembrolizumab.",
    },
    "prostate": {
        "dominant_tme_class": "cold / immunosuppressive",
        "key_immune_features": [
            "Low TMB and minimal immune infiltration",
            "AR signaling suppresses immune response",
            "dMMR/MSI-H rare (~3%) but ICB-responsive subset",
        ],
        "typical_cell_composition": {
            "T_cells_pct": 5, "B_cells_pct": 2, "NK_pct": 2,
            "myeloid_pct": 20, "fibroblast_pct": 35,
        },
        "checkpoint_expression": ["PD-L1 (low)", "CTLA-4 (low)"],
        "treatment_response_pattern": "Generally ICB-resistant; pembrolizumab "
                                      "for dMMR/MSI-H subset. Combination with "
                                      "AR pathway inhibition under study.",
    },
}


# ===================================================================
# PEDIATRIC SINGLE-CELL KNOWLEDGE
# ===================================================================

PEDIATRIC_ALL_IMMUNOPHENOTYPING: Dict[str, Dict[str, Any]] = {
    "pre_B_ALL": {
        "description": "Precursor B-cell ALL (most common pediatric ALL subtype, ~85%)",
        "markers_positive": ["CD19", "CD10", "CD22", "TdT"],
        "markers_negative": ["cytoplasmic mu (IgM)"],
        "scRNA_seq_signature": "PAX5+, EBF1+, CD79A+; arrested at pre-B stage",
    },
    "pro_B_ALL": {
        "description": "Pro-B ALL (early B-lineage, ~5% of B-ALL)",
        "markers_positive": ["CD19", "CD22", "TdT"],
        "markers_negative": ["CD10"],
        "scRNA_seq_signature": "PAX5 low, CD34+; often MLL-rearranged",
    },
    "T_ALL": {
        "description": "T-cell ALL (~15% of pediatric ALL)",
        "markers_positive": ["CD3 (cytoplasmic or surface)", "CD7", "CD5", "TdT"],
        "markers_negative": ["CD19", "CD22 (B-cell markers)"],
        "scRNA_seq_signature": "NOTCH1 pathway activation, TAL1/LMO2 expression",
    },
    "MRD_detection": {
        "description": "Minimal residual disease detection by single-cell approaches",
        "sensitivity": "10^-4 to 10^-5 by flow cytometry and scRNA-seq",
        "clinical_significance": "MRD >0.01% at end-of-induction is adverse prognostic factor",
        "advantage_of_scRNA": "Identifies leukemia-associated immunophenotype shifts that may "
                             "escape standard MFC panels",
    },
}

PEDIATRIC_TUMOR_MICROENVIRONMENT: Dict[str, Dict[str, Any]] = {
    "neuroblastoma": {
        "key_feature": "Schwann cell stroma ratio predicts outcome",
        "tme_characteristics": [
            "Stroma-rich tumors (Schwannian stroma) have favorable prognosis",
            "Stroma-poor tumors are undifferentiated and aggressive",
            "Tumor-associated macrophages correlate with MYCN amplification",
        ],
        "immune_profile": "Variable; MYCN-amplified tumors are immune-cold",
    },
    "medulloblastoma": {
        "key_feature": "Group 3 has immune-cold TME with low T-cell infiltration",
        "tme_by_subgroup": {
            "WNT": "Moderate immune infiltration; best prognosis",
            "SHH": "Variable TME; infant vs adult SHH differ substantially",
            "Group_3": "Immune-cold; minimal CD8+ T cells; worst prognosis",
            "Group_4": "Low-to-moderate immune infiltration",
        },
    },
    "ewing_sarcoma": {
        "key_feature": "Immune desert with minimal CD8+ T cells",
        "tme_characteristics": [
            "Very low TMB (~0.15 mutations/Mb)",
            "Minimal CD8+ T cell infiltration",
            "Predominantly myeloid-derived suppressor cells",
            "EWS-FLI1 fusion may suppress antigen presentation",
        ],
        "immunotherapy_implication": "Poor response to checkpoint inhibitors; "
                                    "CAR-T and bispecific antibodies under investigation",
    },
}

PEDIATRIC_CAR_T_TARGETS: Dict[str, Dict[str, Any]] = {
    "CD19": {
        "target": "CD19",
        "tumor_type": "B-cell ALL",
        "expression_level": ">95% of B-ALL blasts (strong, homogeneous)",
        "fda_approved_products": ["tisagenlecleucel (Kymriah)", "brexucabtagene autoleucel (Tecartus)"],
        "response_rate": "70-90% CR in relapsed/refractory pediatric B-ALL",
        "resistance_mechanism": "CD19 antigen loss in 10-20% of relapses",
    },
    "CD22": {
        "target": "CD22",
        "tumor_type": "B-cell ALL (especially post-CD19 CAR-T relapse)",
        "expression_level": "80-90% of B-ALL blasts (dimmer, heterogeneous → dose-dependent efficacy)",
        "clinical_status": "Multiple clinical trials; dual CD19/CD22 CAR-T in development",
        "response_rate": "60-70% CR; lower durability than CD19 CAR-T alone",
        "resistance_mechanism": "CD22 downregulation and dim expression below CAR threshold",
    },
    "GD2": {
        "target": "GD2 (disialoganglioside)",
        "tumor_type": "Neuroblastoma",
        "expression_level": ">90% of neuroblastoma tumor cells; minimal expression on normal tissues "
                           "(restricted to peripheral nerves, brain)",
        "clinical_status": "Phase I/II trials (e.g., NCT03373097); anti-GD2 antibody (dinutuximab) "
                          "already FDA-approved for high-risk neuroblastoma",
        "challenges": "On-target off-tumor pain (peripheral nerve expression); solid tumor TME barriers",
    },
}
