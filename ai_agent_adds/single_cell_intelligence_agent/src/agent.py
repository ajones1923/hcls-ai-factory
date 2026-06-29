"""Single-Cell Intelligence Agent -- autonomous reasoning across single-cell data.

Implements the plan -> search -> evaluate -> synthesize -> report pattern from the
VAST AI OS AgentEngine model. The agent can:

1. Parse complex multi-part questions about single-cell genomics and spatial biology
2. Plan a search strategy across 12 domain-specific collections
3. Execute multi-collection retrieval via the SingleCellRAGEngine
4. Evaluate evidence quality and completeness
5. Synthesize cross-functional insights with clinical alerts
6. Generate structured reports with single-cell-specific formatting

Mapping to VAST AI OS:
  - AgentEngine entry point: SingleCellAgent.run()
  - Plan -> search_plan()
  - Execute -> rag_engine.query()
  - Reflect -> evaluate_evidence()
  - Report -> generate_report()

Domain coverage:
  - Cell type annotation (reference-based, marker scoring, LLM-assisted consensus)
  - Tumor microenvironment classification (hot/cold/excluded/immunosuppressive)
  - Drug response prediction at cellular resolution
  - Subclonal architecture detection and clonal dynamics
  - Spatial transcriptomics and spatial niche identification
  - Developmental trajectory inference and pseudotime analysis
  - Ligand-receptor interaction mapping and cell communication
  - Cell-type-specific biomarker discovery
  - CAR-T target validation (on-tumor/off-tumor expression profiling)
  - Treatment response monitoring through longitudinal clonal dynamics
  - GPU-accelerated analysis via RAPIDS cuML/cuGraph

Reference databases: Human Cell Atlas, Cell Ontology, CellMarker, PanglaoDB, DepMap

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Dict, List, Optional


# =====================================================================
# ENUMS
# =====================================================================

class SCWorkflowType(str, Enum):
    """Types of single-cell analysis workflows."""
    CELL_TYPE_ANNOTATION = "cell_type_annotation"
    TME_CLASSIFICATION = "tme_classification"
    DRUG_RESPONSE = "drug_response"
    SUBCLONAL_ARCHITECTURE = "subclonal_architecture"
    SPATIAL_ANALYSIS = "spatial_analysis"
    TRAJECTORY_INFERENCE = "trajectory_inference"
    LIGAND_RECEPTOR = "ligand_receptor"
    BIOMARKER_DISCOVERY = "biomarker_discovery"
    CART_TARGET_VALIDATION = "cart_target_validation"
    TREATMENT_MONITORING = "treatment_monitoring"
    GENERAL = "general"


class EvidenceLevel(str, Enum):
    """Clinical evidence hierarchy for single-cell findings."""
    LEVEL_I = "I"        # Large-scale atlas study (HCA, TCGA), prospective clinical
    LEVEL_II = "II"      # Multi-cohort validation, independent replication
    LEVEL_III = "III"    # Single-cohort study, computational prediction validated
    LEVEL_IV = "IV"      # Computational prediction only, case report, expert opinion
    ATLAS = "atlas"      # Human Cell Atlas / major atlas reference data
    BENCHMARK = "benchmark"  # Systematic benchmarking study


class SeverityLevel(str, Enum):
    """Finding severity classification."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INFORMATIONAL = "informational"


class AnalysisModality(str, Enum):
    """Single-cell analysis modalities."""
    SCRNA_SEQ = "scRNA-seq"
    SNRNA_SEQ = "snRNA-seq"
    CITE_SEQ = "CITE-seq"
    SPATIAL_TRANSCRIPTOMICS = "spatial_transcriptomics"
    VISIUM = "Visium"
    MERFISH = "MERFISH"
    SLIDE_SEQ = "Slide-seq"
    MULTIOME = "Multiome"
    SCATAC_SEQ = "scATAC-seq"
    SCRNA_TCR = "scRNA+TCR"
    MASS_CYTOMETRY = "mass_cytometry"


class CellOntologyDomain(str, Enum):
    """High-level Cell Ontology domains."""
    IMMUNE = "immune"
    EPITHELIAL = "epithelial"
    STROMAL = "stromal"
    ENDOTHELIAL = "endothelial"
    NEURONAL = "neuronal"
    STEM_PROGENITOR = "stem_progenitor"
    MALIGNANT = "malignant"


# =====================================================================
# RESPONSE DATACLASS
# =====================================================================

@dataclass
class SCResponse:
    """Complete response from the single-cell intelligence agent.

    Attributes:
        question: Original query text.
        answer: LLM-synthesised answer text.
        results: Ranked search results used for synthesis.
        workflow: SC workflow that was applied.
        confidence: Overall confidence score (0.0 - 1.0).
        citations: Formatted citation list.
        search_time_ms: Total search time in milliseconds.
        collections_searched: Number of collections queried.
        patient_context_used: Whether patient context was injected.
        clinical_alerts: Any critical clinical findings flagged.
        timestamp: ISO 8601 timestamp of response generation.
    """
    question: str = ""
    answer: str = ""
    results: list = field(default_factory=list)
    workflow: Optional[SCWorkflowType] = None
    confidence: float = 0.0
    citations: List[Dict[str, str]] = field(default_factory=list)
    search_time_ms: float = 0.0
    collections_searched: int = 0
    patient_context_used: bool = False
    clinical_alerts: List[str] = field(default_factory=list)
    timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )


# =====================================================================
# SINGLE-CELL SYSTEM PROMPT
# =====================================================================

SC_SYSTEM_PROMPT = """\
You are a single-cell genomics intelligence system within the HCLS AI Factory. \
You analyze scRNA-seq, spatial transcriptomics, and multi-modal single-cell data \
at cellular resolution. You perform cell type annotation using multi-strategy \
consensus (reference-based, marker scoring, LLM-assisted), classify tumor \
microenvironments (hot/cold/excluded/immunosuppressive), predict drug response \
at cellular resolution, detect subclonal architecture, identify spatial niches, \
infer developmental trajectories, map ligand-receptor interactions, discover \
cell-type-specific biomarkers, validate CAR-T targets for on-tumor/off-tumor \
expression, and monitor treatment response through longitudinal clonal dynamics. \
You reference the Human Cell Atlas, Cell Ontology, CellMarker, PanglaoDB, and \
DepMap. You understand GPU-accelerated analysis via RAPIDS cuML/cuGraph. Always \
cite specific cell ontology IDs, marker genes, and evidence levels.

Your responses must adhere to the following standards:

1. **Cell Ontology Citations** -- Always cite cell types using Cell Ontology IDs \
   with clickable links: [CL:0000084](https://www.ebi.ac.uk/ols4/ontologies/cl/classes/CL_0000084) \
   (T cell). When referencing markers, cite the source database: CellMarker, \
   PanglaoDB, or Human Cell Atlas with dataset identifiers. Include canonical \
   marker gene sets and confidence of annotation.

2. **Literature References** -- Cite single-cell studies using PubMed identifiers \
   with clickable links: [PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). \
   Include study type (atlas, cohort, benchmark), sample size (number of cells), \
   and key findings. For clinical trials with scRNA-seq correlatives, cite NCT \
   identifiers: [NCT01234567](https://clinicaltrials.gov/study/NCT01234567).

3. **CRITICAL Findings** -- Flag the following as CRITICAL with prominent visual \
   markers and immediate action recommendations:
   - Aberrant cell populations suggesting malignant transformation
   - On-tumor target expressed on critical normal tissues (CAR-T safety)
   - Immune evasion signatures (PD-L1 high, MHC-I loss, TGF-beta dominance)
   - Drug resistance clones detected at significant frequency (>5%)
   - Spatial exclusion of effector immune cells from tumor core
   - Clonal expansion of T cell clones with exhaustion signatures
   - Loss of tumor suppressor expression in stem/progenitor compartment
   - Therapy-resistant subclones with actionable mutations
   - Neurotoxicity-associated cell states in CAR-T products
   - Unexpected lineage infidelity or trans-differentiation events

4. **Severity Badges** -- Classify all findings using standardised severity levels: \
   [CRITICAL], [HIGH], [MODERATE], [LOW], [INFORMATIONAL]. Place the badge at the \
   start of each finding or recommendation line.

5. **Cell Type Annotation Scoring** -- When performing cell type annotation, show \
   the multi-strategy breakdown:
   - Reference-based: correlation score against HCA/atlas references (0-1.0)
   - Marker scoring: enrichment of canonical markers from CellMarker/PanglaoDB
   - Consensus call: agreement across methods with confidence tier (high/medium/low)
   - Provide Cell Ontology ID, canonical markers, and number of cells in cluster
   - Flag any ambiguous or novel cell states requiring expert review

6. **TME Classification** -- For tumor microenvironment queries, provide structured \
   classification using the four-category framework:
   - Hot (immune-inflamed): CD8+ T cell infiltration, IFN-gamma signature, PD-L1+
   - Cold (immune-desert): absence of T cell infiltration, low MHC-I
   - Excluded: immune cells present at tumor margin but not core, TGF-beta signature
   - Immunosuppressive: Treg/MDSC enrichment, IL-10/TGF-beta high, M2 macrophages
   Include spatial context when available and treatment implications.

7. **Structured Formatting** -- Organise responses with clear sections: \
   Cell Population Summary, Marker Gene Analysis, Spatial Context, \
   TME Classification, Drug Response Predictions, Trajectory Analysis, \
   Ligand-Receptor Interactions, Clinical Implications, and \
   Methodology Notes. Use bullet points and numbered lists for \
   actionable items.

8. **Genomic Cross-Reference** -- When relevant mutations or CNVs are detected \
   in single-cell data, cross-reference with the genomic_evidence collection. \
   Discuss clonal architecture (founder vs. subclonal), allele frequencies, \
   therapy implications, and integration with bulk WGS/WES findings. Reference \
   DepMap dependency scores for druggability assessment.

9. **Method Recommendations** -- For analytical questions, recommend specific \
   computational methods with evidence:
   - Clustering: Leiden (preferred for large datasets) vs. Louvain vs. spectral
   - Integration: Harmony, scVI, SCVI, scANVI for batch correction
   - Trajectory: Monocle3, scVelo (RNA velocity), PAGA, CellRank
   - Spatial: SpatialDE, SPARK, cell2location, Tangram, CellCharter
   - Communication: CellChat, LIANA+, NicheNet, Commot
   - Annotation: CellTypist, Azimuth, scArches, SingleR
   - GPU acceleration: RAPIDS cuML/cuGraph for scalability
   Always cite benchmarking studies supporting method selection.

10. **Limitations** -- You are a single-cell genomics intelligence tool. You \
    do NOT replace expert bioinformaticians, pathologists, or oncologists. \
    Cell type annotations are computational predictions requiring experimental \
    validation. Drug response predictions from scRNA-seq are hypothesis-generating \
    and require functional validation. Spatial inferences depend on capture \
    technology resolution. Explicitly state when evidence is limited, when \
    experimental validation is needed, or when specialist consultation \
    (computational biology, pathology, immunology) is recommended."""


# =====================================================================
# WORKFLOW-SPECIFIC COLLECTION BOOST WEIGHTS
# =====================================================================
# Maps each SCWorkflowType to collection weight overrides (multipliers).
# Collections not listed retain their base weight (1.0x). Values > 1.0
# boost the collection; values < 1.0 would suppress it.

WORKFLOW_COLLECTION_BOOST: Dict[SCWorkflowType, Dict[str, float]] = {

    # -- Cell Type Annotation -------------------------------------------
    SCWorkflowType.CELL_TYPE_ANNOTATION: {
        "sc_cell_types": 2.5,
        "sc_markers": 2.2,
        "sc_literature": 1.5,
        "sc_methods": 1.3,
        "sc_datasets": 1.2,
        "sc_trajectories": 0.8,
        "sc_spatial": 0.8,
        "sc_tme": 0.7,
        "sc_drug_response": 0.5,
        "sc_pathways": 0.8,
        "sc_clinical": 0.6,
        "genomic_evidence": 0.7,
    },

    # -- TME Classification ---------------------------------------------
    SCWorkflowType.TME_CLASSIFICATION: {
        "sc_tme": 2.5,
        "sc_cell_types": 2.0,
        "sc_spatial": 1.8,
        "sc_markers": 1.5,
        "sc_literature": 1.5,
        "sc_drug_response": 1.3,
        "sc_clinical": 1.2,
        "sc_pathways": 1.0,
        "sc_methods": 0.8,
        "sc_datasets": 0.7,
        "sc_trajectories": 0.8,
        "genomic_evidence": 0.8,
    },

    # -- Drug Response --------------------------------------------------
    SCWorkflowType.DRUG_RESPONSE: {
        "sc_drug_response": 2.5,
        "sc_clinical": 2.0,
        "sc_tme": 1.5,
        "sc_cell_types": 1.3,
        "sc_pathways": 1.5,
        "sc_literature": 1.3,
        "sc_markers": 1.0,
        "sc_methods": 0.8,
        "sc_trajectories": 0.8,
        "sc_spatial": 0.8,
        "sc_datasets": 0.7,
        "genomic_evidence": 1.2,
    },

    # -- Subclonal Architecture -----------------------------------------
    SCWorkflowType.SUBCLONAL_ARCHITECTURE: {
        "sc_cell_types": 2.0,
        "sc_trajectories": 2.0,
        "sc_clinical": 1.5,
        "sc_literature": 1.5,
        "sc_pathways": 1.3,
        "sc_drug_response": 1.2,
        "sc_tme": 1.0,
        "sc_markers": 1.0,
        "sc_methods": 1.2,
        "sc_spatial": 0.8,
        "sc_datasets": 0.8,
        "genomic_evidence": 2.0,
    },

    # -- Spatial Analysis -----------------------------------------------
    SCWorkflowType.SPATIAL_ANALYSIS: {
        "sc_spatial": 2.5,
        "sc_tme": 2.0,
        "sc_cell_types": 1.8,
        "sc_markers": 1.3,
        "sc_methods": 1.5,
        "sc_literature": 1.3,
        "sc_datasets": 1.2,
        "sc_pathways": 1.0,
        "sc_drug_response": 0.8,
        "sc_trajectories": 0.8,
        "sc_clinical": 0.7,
        "genomic_evidence": 0.7,
    },

    # -- Trajectory Inference -------------------------------------------
    SCWorkflowType.TRAJECTORY_INFERENCE: {
        "sc_trajectories": 2.5,
        "sc_cell_types": 2.0,
        "sc_markers": 1.5,
        "sc_methods": 1.8,
        "sc_literature": 1.3,
        "sc_pathways": 1.2,
        "sc_datasets": 1.0,
        "sc_spatial": 0.8,
        "sc_tme": 0.7,
        "sc_drug_response": 0.7,
        "sc_clinical": 0.6,
        "genomic_evidence": 0.8,
    },

    # -- Ligand-Receptor ------------------------------------------------
    SCWorkflowType.LIGAND_RECEPTOR: {
        "sc_pathways": 2.5,
        "sc_cell_types": 2.0,
        "sc_spatial": 1.8,
        "sc_tme": 1.5,
        "sc_markers": 1.3,
        "sc_literature": 1.3,
        "sc_methods": 1.2,
        "sc_drug_response": 1.0,
        "sc_datasets": 0.8,
        "sc_trajectories": 0.8,
        "sc_clinical": 0.7,
        "genomic_evidence": 0.7,
    },

    # -- Biomarker Discovery --------------------------------------------
    SCWorkflowType.BIOMARKER_DISCOVERY: {
        "sc_markers": 2.5,
        "sc_cell_types": 2.0,
        "sc_clinical": 1.8,
        "sc_literature": 1.5,
        "sc_drug_response": 1.3,
        "sc_tme": 1.2,
        "sc_pathways": 1.2,
        "sc_datasets": 1.0,
        "sc_methods": 1.0,
        "sc_spatial": 0.8,
        "sc_trajectories": 0.8,
        "genomic_evidence": 1.2,
    },

    # -- CAR-T Target Validation ----------------------------------------
    SCWorkflowType.CART_TARGET_VALIDATION: {
        "sc_cell_types": 2.5,
        "sc_markers": 2.2,
        "sc_tme": 2.0,
        "sc_clinical": 1.8,
        "sc_drug_response": 1.5,
        "sc_spatial": 1.3,
        "sc_literature": 1.3,
        "sc_pathways": 1.0,
        "sc_datasets": 1.0,
        "sc_methods": 0.8,
        "sc_trajectories": 0.8,
        "genomic_evidence": 1.5,
    },

    # -- Treatment Monitoring -------------------------------------------
    SCWorkflowType.TREATMENT_MONITORING: {
        "sc_clinical": 2.5,
        "sc_cell_types": 2.0,
        "sc_tme": 1.8,
        "sc_drug_response": 1.8,
        "sc_trajectories": 1.5,
        "sc_markers": 1.3,
        "sc_literature": 1.3,
        "sc_pathways": 1.0,
        "sc_spatial": 1.0,
        "sc_methods": 0.8,
        "sc_datasets": 0.7,
        "genomic_evidence": 1.2,
    },

    # -- General (balanced across all collections) ----------------------
    SCWorkflowType.GENERAL: {
        "sc_cell_types": 1.2,
        "sc_markers": 1.1,
        "sc_literature": 1.2,
        "sc_methods": 1.0,
        "sc_datasets": 1.0,
        "sc_spatial": 1.0,
        "sc_tme": 1.0,
        "sc_drug_response": 0.9,
        "sc_trajectories": 1.0,
        "sc_pathways": 1.0,
        "sc_clinical": 0.9,
        "genomic_evidence": 0.8,
    },
}


# =====================================================================
# KNOWLEDGE DOMAIN DICTIONARIES
# =====================================================================
# Comprehensive single-cell knowledge for entity detection and context
# enrichment. Used by the agent's search_plan() to identify entities
# in user queries and map them to workflows.

SC_CONDITIONS: Dict[str, Dict[str, object]] = {

    # -- Solid Tumors (scRNA-seq context) --------------------------------
    "non-small cell lung cancer": {
        "aliases": ["nsclc", "lung adenocarcinoma", "luad", "lung squamous",
                    "lusc", "lung cancer"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.DRUG_RESPONSE],
        "search_terms": ["EGFR", "ALK", "KRAS G12C", "PD-L1", "TPS",
                        "immune checkpoint", "tumor-infiltrating lymphocytes",
                        "alveolar type 2", "club cell"],
    },
    "breast cancer": {
        "aliases": ["brca", "triple negative breast cancer", "tnbc",
                    "er positive", "her2 positive", "luminal a", "luminal b"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.SUBCLONAL_ARCHITECTURE],
        "search_terms": ["ER", "PR", "HER2", "BRCA1", "BRCA2", "luminal",
                        "basal-like", "tumor heterogeneity", "cancer-associated fibroblasts"],
    },
    "colorectal cancer": {
        "aliases": ["crc", "colon cancer", "rectal cancer", "colorectal adenocarcinoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.DRUG_RESPONSE],
        "search_terms": ["MSI-H", "microsatellite instability", "APC", "KRAS",
                        "WNT pathway", "stem cell niche", "goblet cell",
                        "enterocyte", "consensus molecular subtype"],
    },
    "melanoma": {
        "aliases": ["cutaneous melanoma", "metastatic melanoma", "uveal melanoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.TREATMENT_MONITORING],
        "search_terms": ["BRAF V600E", "NRAS", "immune checkpoint inhibitor",
                        "melanocyte", "neural crest", "T cell exhaustion",
                        "immunotherapy response", "MHC-I expression"],
    },
    "glioblastoma": {
        "aliases": ["gbm", "glioblastoma multiforme", "high-grade glioma",
                    "idh-wildtype glioblastoma"],
        "workflows": [SCWorkflowType.SUBCLONAL_ARCHITECTURE, SCWorkflowType.TME_CLASSIFICATION],
        "search_terms": ["glioma stem cell", "mesenchymal", "proneural",
                        "classical", "neural", "MGMT methylation",
                        "IDH mutation", "cell state plasticity",
                        "tumor-associated macrophage"],
    },
    "pancreatic ductal adenocarcinoma": {
        "aliases": ["pdac", "pancreatic cancer", "pancreatic adenocarcinoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.SPATIAL_ANALYSIS],
        "search_terms": ["KRAS", "stellate cell", "desmoplastic stroma",
                        "immune exclusion", "ductal cell", "acinar cell",
                        "cancer-associated fibroblast", "myofibroblast"],
    },
    "hepatocellular carcinoma": {
        "aliases": ["hcc", "liver cancer", "hepatocellular"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.BIOMARKER_DISCOVERY],
        "search_terms": ["AFP", "hepatocyte", "Kupffer cell", "sinusoidal",
                        "cholangiocyte", "immune evasion", "HBV", "HCV",
                        "tumor-associated macrophage"],
    },
    "renal cell carcinoma": {
        "aliases": ["rcc", "clear cell rcc", "kidney cancer", "ccrcc"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.DRUG_RESPONSE],
        "search_terms": ["VHL", "HIF pathway", "sunitinib", "nivolumab",
                        "proximal tubule", "immune infiltrate",
                        "angiogenesis", "PBRM1"],
    },
    "ovarian cancer": {
        "aliases": ["ovarian carcinoma", "high-grade serous ovarian",
                    "hgsoc", "ovarian serous"],
        "workflows": [SCWorkflowType.SUBCLONAL_ARCHITECTURE, SCWorkflowType.DRUG_RESPONSE],
        "search_terms": ["BRCA1", "BRCA2", "TP53", "homologous recombination",
                        "PARP inhibitor", "platinum resistance",
                        "fallopian tube epithelium", "mesothelial"],
    },
    "head and neck squamous cell carcinoma": {
        "aliases": ["hnscc", "head and neck cancer", "oral squamous cell carcinoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.SPATIAL_ANALYSIS],
        "search_terms": ["HPV", "PD-L1", "keratinocyte", "basal cell",
                        "T cell exhaustion", "spatial immune context",
                        "cetuximab", "pembrolizumab"],
    },

    # -- Hematologic Malignancies ----------------------------------------
    "acute myeloid leukemia": {
        "aliases": ["aml", "acute myeloid leukaemia"],
        "workflows": [SCWorkflowType.SUBCLONAL_ARCHITECTURE, SCWorkflowType.TRAJECTORY_INFERENCE],
        "search_terms": ["leukemia stem cell", "FLT3-ITD", "NPM1",
                        "differentiation hierarchy", "blast cell",
                        "clonal evolution", "MRD", "HSC"],
    },
    "diffuse large b-cell lymphoma": {
        "aliases": ["dlbcl", "large b cell lymphoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.CART_TARGET_VALIDATION],
        "search_terms": ["CD19", "CD20", "GCB", "ABC", "CAR-T",
                        "cell of origin", "double hit", "MYC", "BCL2"],
    },
    "multiple myeloma": {
        "aliases": ["myeloma", "mm", "plasma cell myeloma"],
        "workflows": [SCWorkflowType.SUBCLONAL_ARCHITECTURE, SCWorkflowType.DRUG_RESPONSE],
        "search_terms": ["BCMA", "plasma cell", "bone marrow microenvironment",
                        "clonal evolution", "GPRC5D", "FcRH5",
                        "immunoglobulin", "t(4;14)", "del(17p)"],
    },
    "b-cell acute lymphoblastic leukemia": {
        "aliases": ["b-all", "b-cell all", "acute lymphoblastic leukemia",
                    "all"],
        "workflows": [SCWorkflowType.CART_TARGET_VALIDATION, SCWorkflowType.TREATMENT_MONITORING],
        "search_terms": ["CD19", "CD22", "BLIN-CAR", "MRD",
                        "pre-B cell", "Philadelphia chromosome",
                        "BCR-ABL", "antigen escape"],
    },

    # -- Autoimmune Diseases --------------------------------------------
    "systemic lupus erythematosus": {
        "aliases": ["sle", "lupus"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.BIOMARKER_DISCOVERY],
        "search_terms": ["interferon signature", "plasmablast", "B cell",
                        "T follicular helper", "renal", "nephritis",
                        "autoantibody", "type I interferon"],
    },
    "rheumatoid arthritis": {
        "aliases": ["ra", "rheumatoid"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.LIGAND_RECEPTOR],
        "search_terms": ["synovial fibroblast", "macrophage", "T cell",
                        "B cell", "TNF", "IL-6", "RANKL",
                        "synovium", "pannus"],
    },
    "inflammatory bowel disease": {
        "aliases": ["ibd", "crohn's disease", "crohns", "ulcerative colitis",
                    "uc"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.SPATIAL_ANALYSIS],
        "search_terms": ["intestinal epithelium", "goblet cell", "Paneth cell",
                        "lamina propria", "Th17", "ILC3", "TNF",
                        "IL-23", "anti-TNF response"],
    },
    "multiple sclerosis": {
        "aliases": ["ms", "relapsing-remitting ms", "rrms"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.BIOMARKER_DISCOVERY],
        "search_terms": ["oligodendrocyte", "microglia", "astrocyte",
                        "demyelination", "CSF", "B cell",
                        "CD20", "ocrelizumab"],
    },

    # -- Developmental / Normal Biology ---------------------------------
    "normal hematopoiesis": {
        "aliases": ["hematopoiesis", "hsc differentiation", "bone marrow",
                    "myelopoiesis", "lymphopoiesis"],
        "workflows": [SCWorkflowType.TRAJECTORY_INFERENCE, SCWorkflowType.CELL_TYPE_ANNOTATION],
        "search_terms": ["HSC", "MPP", "CMP", "GMP", "MEP",
                        "LMPP", "CLP", "erythropoiesis",
                        "granulopoiesis", "differentiation"],
    },
    "embryonic development": {
        "aliases": ["embryogenesis", "organogenesis", "gastrulation",
                    "fetal development"],
        "workflows": [SCWorkflowType.TRAJECTORY_INFERENCE, SCWorkflowType.SPATIAL_ANALYSIS],
        "search_terms": ["germ layer", "mesoderm", "endoderm", "ectoderm",
                        "neural crest", "somite", "lineage specification",
                        "Wnt", "BMP", "Notch"],
    },
    "neurogenesis": {
        "aliases": ["neural development", "cortical development",
                    "brain organoid", "neural progenitor"],
        "workflows": [SCWorkflowType.TRAJECTORY_INFERENCE, SCWorkflowType.CELL_TYPE_ANNOTATION],
        "search_terms": ["radial glia", "intermediate progenitor",
                        "cortical neuron", "interneuron", "oligodendrocyte precursor",
                        "astrocyte", "neural stem cell", "Notch", "SHH"],
    },

    # -- Fibrotic / Degenerative ----------------------------------------
    "idiopathic pulmonary fibrosis": {
        "aliases": ["ipf", "pulmonary fibrosis", "lung fibrosis"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.LIGAND_RECEPTOR],
        "search_terms": ["alveolar type 2", "myofibroblast", "basal cell",
                        "aberrant basaloid", "TGF-beta", "WNT",
                        "honeycombing", "fibroblast focus"],
    },
    "alzheimer disease": {
        "aliases": ["alzheimers", "alzheimer's disease", "ad"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.SPATIAL_ANALYSIS],
        "search_terms": ["microglia", "DAM", "disease-associated microglia",
                        "astrocyte", "oligodendrocyte", "neuron loss",
                        "amyloid", "tau", "TREM2", "APOE"],
    },

    # -- COVID-19 / Infectious Disease ----------------------------------
    "covid-19": {
        "aliases": ["sars-cov-2", "coronavirus", "covid"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.BIOMARKER_DISCOVERY],
        "search_terms": ["ACE2", "TMPRSS2", "alveolar macrophage",
                        "cytokine storm", "monocyte", "neutrophil",
                        "T cell exhaustion", "interferon response",
                        "BALF", "nasopharyngeal"],
    },

    # -- Additional Solid Tumors -------------------------------------------
    "bladder cancer": {
        "aliases": ["urothelial carcinoma", "bladder urothelial", "mibc",
                    "nmibc", "transitional cell carcinoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.DRUG_RESPONSE],
        "search_terms": ["FGFR3", "urothelium", "basal", "luminal",
                        "neuroendocrine", "immune checkpoint",
                        "BCG response", "umbrella cell"],
    },
    "prostate cancer": {
        "aliases": ["prostate adenocarcinoma", "crpc", "castration-resistant prostate",
                    "prostate carcinoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.SUBCLONAL_ARCHITECTURE],
        "search_terms": ["AR", "androgen receptor", "PSA", "KLK3",
                        "luminal", "basal", "neuroendocrine differentiation",
                        "TMPRSS2-ERG", "PTEN loss"],
    },
    "thyroid cancer": {
        "aliases": ["thyroid carcinoma", "papillary thyroid cancer", "ptc",
                    "anaplastic thyroid cancer"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.TME_CLASSIFICATION],
        "search_terms": ["BRAF V600E", "RAS", "thyrocyte", "follicular cell",
                        "RET/PTC", "PAX8", "NIS", "thyroglobulin",
                        "dedifferentiation"],
    },
    "endometrial cancer": {
        "aliases": ["endometrial carcinoma", "uterine cancer",
                    "endometrial adenocarcinoma", "uterine carcinoma"],
        "workflows": [SCWorkflowType.TME_CLASSIFICATION, SCWorkflowType.BIOMARKER_DISCOVERY],
        "search_terms": ["POLE", "MSI-H", "p53", "PTEN", "endometrial epithelium",
                        "stromal cell", "hormone receptor", "immune infiltrate",
                        "copy number high"],
    },

    # -- Additional Hematologic Malignancies --------------------------------
    "chronic lymphocytic leukemia": {
        "aliases": ["cll", "chronic lymphocytic leukaemia", "sll",
                    "small lymphocytic lymphoma"],
        "workflows": [SCWorkflowType.SUBCLONAL_ARCHITECTURE, SCWorkflowType.TREATMENT_MONITORING],
        "search_terms": ["CD5", "CD23", "IGHV mutation", "del(13q)",
                        "del(17p)", "TP53", "BTK inhibitor", "ibrutinib",
                        "venetoclax", "CLL cell evolution"],
    },
    "myelodysplastic syndrome": {
        "aliases": ["mds", "myelodysplasia", "refractory anemia",
                    "myelodysplastic neoplasm"],
        "workflows": [SCWorkflowType.TRAJECTORY_INFERENCE, SCWorkflowType.CELL_TYPE_ANNOTATION],
        "search_terms": ["SF3B1", "TET2", "DNMT3A", "ASXL1",
                        "ineffective hematopoiesis", "dysplastic lineage",
                        "blast proportion", "ring sideroblast",
                        "clonal hematopoiesis"],
    },

    # -- Autoimmune / Inflammatory -----------------------------------------
    "type 1 diabetes": {
        "aliases": ["t1d", "type1 diabetes", "insulin-dependent diabetes",
                    "iddm", "autoimmune diabetes"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.BIOMARKER_DISCOVERY],
        "search_terms": ["beta cell", "islet cell", "insulitis",
                        "autoantibody", "GAD65", "IA-2", "CD8 autoreactive",
                        "islet destruction", "pancreatic islet"],
    },
    "asthma": {
        "aliases": ["allergic asthma", "eosinophilic asthma", "severe asthma",
                    "airway hyperresponsiveness"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.LIGAND_RECEPTOR],
        "search_terms": ["airway remodeling", "goblet cell metaplasia",
                        "eosinophil", "ILC2", "IL-4", "IL-5", "IL-13",
                        "mast cell", "airway smooth muscle",
                        "bronchial epithelium"],
    },

    # -- Fibrotic / Degenerative / Regenerative ----------------------------
    "liver fibrosis": {
        "aliases": ["hepatic fibrosis", "liver cirrhosis", "nash fibrosis",
                    "hepatic stellate cell activation"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.LIGAND_RECEPTOR],
        "search_terms": ["stellate cell activation", "myofibroblast",
                        "collagen deposition", "TGF-beta", "PDGFR",
                        "hepatocyte injury", "Kupffer cell",
                        "portal fibroblast", "ACTA2"],
    },
    "kidney disease": {
        "aliases": ["chronic kidney disease", "ckd", "glomerulonephritis",
                    "diabetic nephropathy", "renal fibrosis"],
        "workflows": [SCWorkflowType.CELL_TYPE_ANNOTATION, SCWorkflowType.SPATIAL_ANALYSIS],
        "search_terms": ["podocyte loss", "tubular injury", "NPHS1",
                        "NPHS2", "proximal tubule", "collecting duct",
                        "mesangial cell", "glomerular endothelial",
                        "interstitial fibrosis"],
    },
    "cardiac regeneration": {
        "aliases": ["heart regeneration", "cardiac repair",
                    "myocardial regeneration", "cardiomyocyte renewal"],
        "workflows": [SCWorkflowType.TRAJECTORY_INFERENCE, SCWorkflowType.CELL_TYPE_ANNOTATION],
        "search_terms": ["cardiomyocyte proliferation", "cardiac progenitor",
                        "TNNT2", "NKX2-5", "cell cycle re-entry",
                        "fibrotic scar", "cardiac fibroblast",
                        "endocardial cell", "epicardial cell"],
    },
    "wound healing": {
        "aliases": ["tissue repair", "skin wound healing", "fibroblast dynamics",
                    "wound repair", "tissue regeneration"],
        "workflows": [SCWorkflowType.TRAJECTORY_INFERENCE, SCWorkflowType.LIGAND_RECEPTOR],
        "search_terms": ["fibroblast activation", "myofibroblast differentiation",
                        "keratinocyte migration", "re-epithelialization",
                        "granulation tissue", "macrophage polarization",
                        "TGF-beta", "collagen remodeling", "angiogenesis"],
    },
}


SC_CELL_TYPES: Dict[str, Dict[str, object]] = {

    "cd8_t_cell": {
        "full_name": "CD8+ Cytotoxic T Lymphocyte",
        "cell_ontology_id": "CL:0000625",
        "canonical_markers": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "NKG7"],
        "subtypes": ["naive", "effector", "memory", "exhausted", "tissue-resident"],
        "domain": "immune",
    },
    "cd4_t_cell": {
        "full_name": "CD4+ T Helper Cell",
        "cell_ontology_id": "CL:0000624",
        "canonical_markers": ["CD4", "IL7R", "TCF7", "CCR7", "FOXP3", "IL2RA"],
        "subtypes": ["naive", "Th1", "Th2", "Th17", "Tfh", "Treg"],
        "domain": "immune",
    },
    "regulatory_t_cell": {
        "full_name": "Regulatory T Cell (Treg)",
        "cell_ontology_id": "CL:0000815",
        "canonical_markers": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "CD4"],
        "subtypes": ["thymic", "peripheral", "effector Treg"],
        "domain": "immune",
    },
    "b_cell": {
        "full_name": "B Lymphocyte",
        "cell_ontology_id": "CL:0000236",
        "canonical_markers": ["CD19", "CD79A", "MS4A1", "PAX5", "CD20"],
        "subtypes": ["naive", "memory", "germinal center", "plasmablast", "plasma cell"],
        "domain": "immune",
    },
    "nk_cell": {
        "full_name": "Natural Killer Cell",
        "cell_ontology_id": "CL:0000623",
        "canonical_markers": ["NKG7", "GNLY", "KLRD1", "NCAM1", "NCR1", "KLRB1"],
        "subtypes": ["CD56bright", "CD56dim", "adaptive NK"],
        "domain": "immune",
    },
    "monocyte": {
        "full_name": "Monocyte",
        "cell_ontology_id": "CL:0000576",
        "canonical_markers": ["CD14", "LYZ", "S100A8", "S100A9", "FCGR3A", "VCAN"],
        "subtypes": ["classical CD14+", "intermediate", "non-classical CD16+"],
        "domain": "immune",
    },
    "macrophage": {
        "full_name": "Macrophage",
        "cell_ontology_id": "CL:0000235",
        "canonical_markers": ["CD68", "CD163", "MARCO", "MSR1", "MRC1", "CSF1R"],
        "subtypes": ["M1", "M2", "tumor-associated", "tissue-resident", "alveolar"],
        "domain": "immune",
    },
    "dendritic_cell": {
        "full_name": "Dendritic Cell",
        "cell_ontology_id": "CL:0000451",
        "canonical_markers": ["ITGAX", "HLA-DRA", "FLT3", "CLEC9A", "CD1C", "LAMP3"],
        "subtypes": ["cDC1", "cDC2", "pDC", "mature DC", "migratory DC"],
        "domain": "immune",
    },
    "neutrophil": {
        "full_name": "Neutrophil",
        "cell_ontology_id": "CL:0000775",
        "canonical_markers": ["CSF3R", "FCGR3B", "CXCR2", "S100A8", "S100A9", "MMP9"],
        "subtypes": ["mature", "immature", "low-density", "tumor-associated"],
        "domain": "immune",
    },
    "mast_cell": {
        "full_name": "Mast Cell",
        "cell_ontology_id": "CL:0000097",
        "canonical_markers": ["KIT", "TPSAB1", "TPSB2", "CPA3", "HPGDS", "HDC"],
        "subtypes": ["mucosal", "connective tissue"],
        "domain": "immune",
    },
    "fibroblast": {
        "full_name": "Fibroblast",
        "cell_ontology_id": "CL:0000057",
        "canonical_markers": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "PDGFRA"],
        "subtypes": ["myofibroblast", "inflammatory CAF", "myCAF", "iCAF", "apCAF"],
        "domain": "stromal",
    },
    "endothelial_cell": {
        "full_name": "Endothelial Cell",
        "cell_ontology_id": "CL:0000115",
        "canonical_markers": ["PECAM1", "VWF", "CDH5", "ERG", "FLT1", "KDR"],
        "subtypes": ["arterial", "venous", "capillary", "lymphatic", "tip cell"],
        "domain": "endothelial",
    },
    "epithelial_cell": {
        "full_name": "Epithelial Cell",
        "cell_ontology_id": "CL:0000066",
        "canonical_markers": ["EPCAM", "KRT18", "KRT19", "CDH1", "KRT8"],
        "subtypes": ["basal", "luminal", "secretory", "ciliated", "AT2", "AT1"],
        "domain": "epithelial",
    },
    "hematopoietic_stem_cell": {
        "full_name": "Hematopoietic Stem Cell",
        "cell_ontology_id": "CL:0000037",
        "canonical_markers": ["CD34", "KIT", "THY1", "CRHBP", "HLF", "AVP"],
        "subtypes": ["long-term HSC", "short-term HSC", "MPP"],
        "domain": "stem_progenitor",
    },
    "cancer_stem_cell": {
        "full_name": "Cancer Stem Cell / Tumor-Initiating Cell",
        "cell_ontology_id": "CL:0001064",
        "canonical_markers": ["CD44", "ALDH1A1", "PROM1", "SOX2", "NANOG", "OCT4"],
        "subtypes": ["quiescent", "proliferative", "drug-resistant"],
        "domain": "malignant",
    },
    "oligodendrocyte": {
        "full_name": "Oligodendrocyte",
        "cell_ontology_id": "CL:0000128",
        "canonical_markers": ["MBP", "MOG", "PLP1", "MAG", "OLIG2", "SOX10"],
        "subtypes": ["OPC", "pre-myelinating", "mature myelinating"],
        "domain": "neuronal",
    },
    "microglia": {
        "full_name": "Microglia",
        "cell_ontology_id": "CL:0000129",
        "canonical_markers": ["TMEM119", "P2RY12", "CX3CR1", "AIF1", "CSF1R", "TREM2"],
        "subtypes": ["homeostatic", "DAM", "activated", "phagocytic"],
        "domain": "neuronal",
    },
    "astrocyte": {
        "full_name": "Astrocyte",
        "cell_ontology_id": "CL:0000127",
        "canonical_markers": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1", "S100B", "GJA1"],
        "subtypes": ["protoplasmic", "fibrous", "reactive", "disease-associated"],
        "domain": "neuronal",
    },
    "hepatocyte": {
        "full_name": "Hepatocyte",
        "cell_ontology_id": "CL:0000182",
        "canonical_markers": ["ALB", "APOB", "HNF4A", "CYP3A4", "TTR", "SERPINA1"],
        "subtypes": ["periportal", "pericentral", "midzonal"],
        "domain": "epithelial",
    },
    "cardiomyocyte": {
        "full_name": "Cardiomyocyte",
        "cell_ontology_id": "CL:0000746",
        "canonical_markers": ["TNNT2", "MYH7", "MYL2", "ACTC1", "RYR2", "TTN"],
        "subtypes": ["atrial", "ventricular", "nodal"],
        "domain": "stromal",
    },
    "alveolar_type2": {
        "full_name": "Alveolar Type 2 Cell",
        "cell_ontology_id": "CL:0002063",
        "canonical_markers": ["SFTPC", "SFTPB", "ABCA3", "NKX2-1", "LAMP3", "ETV5"],
        "subtypes": ["normal AT2", "transitional", "aberrant basaloid"],
        "domain": "epithelial",
    },

    # -- Dendritic cell subtypes -------------------------------------------
    "plasmacytoid_dc": {
        "full_name": "Plasmacytoid Dendritic Cell",
        "cell_ontology_id": "CL:0000784",
        "canonical_markers": ["CLEC4C", "IL3RA", "IRF7", "TCF4", "LILRA4", "JCHAIN"],
        "subtypes": ["resting pDC", "activated pDC", "IFN-producing pDC"],
        "domain": "immune",
    },
    "cdc1": {
        "full_name": "Conventional Dendritic Cell Type 1",
        "cell_ontology_id": "CL:0002394",
        "canonical_markers": ["CLEC9A", "XCR1", "BATF3", "IRF8", "IDO1", "CADM1"],
        "subtypes": ["immature cDC1", "mature cDC1", "migratory cDC1"],
        "domain": "immune",
    },

    # -- Innate lymphoid cells ---------------------------------------------
    "ilc1": {
        "full_name": "Innate Lymphoid Cell Type 1",
        "cell_ontology_id": "CL:0001077",
        "canonical_markers": ["TBX21", "IFNG", "IL12RB2", "NCR1", "KLRB1", "IL7R"],
        "subtypes": ["tissue-resident ILC1", "inflammatory ILC1"],
        "domain": "immune",
    },
    "ilc2": {
        "full_name": "Innate Lymphoid Cell Type 2",
        "cell_ontology_id": "CL:0001069",
        "canonical_markers": ["GATA3", "IL5", "IL13", "PTGDR2", "IL1RL1", "KLRG1"],
        "subtypes": ["natural ILC2", "inflammatory ILC2", "tissue-resident ILC2"],
        "domain": "immune",
    },
    "ilc3": {
        "full_name": "Innate Lymphoid Cell Type 3",
        "cell_ontology_id": "CL:0001071",
        "canonical_markers": ["RORC", "IL17A", "IL22", "IL23R", "NCR2", "KIT"],
        "subtypes": ["NCR+ ILC3", "NCR- ILC3", "LTi-like ILC3"],
        "domain": "immune",
    },

    # -- Unconventional T cells --------------------------------------------
    "gamma_delta_t_cell": {
        "full_name": "Gamma-Delta T Cell",
        "cell_ontology_id": "CL:0000798",
        "canonical_markers": ["TRGV9", "TRDV2", "CD3E", "TRDC", "KLRC1", "NKG7"],
        "subtypes": ["Vgamma9Vdelta2", "Vdelta1", "tissue-resident gdT"],
        "domain": "immune",
    },
    "mait_cell": {
        "full_name": "Mucosal-Associated Invariant T Cell",
        "cell_ontology_id": "CL:0000940",
        "canonical_markers": ["TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16", "IL18R1", "DPP4"],
        "subtypes": ["MAIT1", "MAIT17", "tissue-resident MAIT"],
        "domain": "immune",
    },

    # -- Stromal and structural cells --------------------------------------
    "mast_cell_v2": {
        "full_name": "Mast Cell (Extended Markers)",
        "cell_ontology_id": "CL:0000097",
        "canonical_markers": ["KIT", "FCER1A", "TPSB2", "CPA3", "HPGDS", "HDC"],
        "subtypes": ["MCT (tryptase+)", "MCTC (tryptase+chymase+)", "mucosal", "connective tissue"],
        "domain": "immune",
    },
    "pericyte": {
        "full_name": "Pericyte",
        "cell_ontology_id": "CL:0000669",
        "canonical_markers": ["RGS5", "PDGFRB", "NOTCH3", "ACTA2", "DES", "CSPG4"],
        "subtypes": ["arteriolar pericyte", "capillary pericyte", "venular pericyte"],
        "domain": "stromal",
    },
    "cancer_associated_fibroblast": {
        "full_name": "Cancer-Associated Fibroblast",
        "cell_ontology_id": "CL:0000057",
        "canonical_markers": ["FAP", "POSTN", "COL1A1", "COL3A1", "ACTA2", "PDPN"],
        "subtypes": ["myCAF", "iCAF", "apCAF", "desmoplastic CAF", "matrix CAF"],
        "domain": "stromal",
    },
}


SC_BIOMARKERS: Dict[str, Dict[str, str]] = {
    "cd19-expression": {
        "full_name": "CD19 Surface Expression (scRNA-seq / CITE-seq)",
        "assay": "scRNA-seq, CITE-seq, flow cytometry",
        "significance": "Primary CAR-T target in B-ALL and DLBCL; scRNA-seq enables "
                        "detection of antigen-low/negative subclones that may escape "
                        "CAR-T therapy; CITE-seq validates protein-level expression; "
                        "critical for predicting CD19-negative relapse risk",
        "workflows": "cart_target_validation,treatment_monitoring",
    },
    "bcma-expression": {
        "full_name": "BCMA (TNFRSF17) Expression",
        "assay": "scRNA-seq, CITE-seq, flow cytometry",
        "significance": "CAR-T and bispecific antibody target in multiple myeloma; "
                        "single-cell quantification reveals expression heterogeneity "
                        "across plasma cell subclones; soluble BCMA shedding detectable "
                        "at cellular level",
        "workflows": "cart_target_validation,drug_response",
    },
    "pd-l1-expression": {
        "full_name": "PD-L1 (CD274) Expression by Cell Type",
        "assay": "scRNA-seq, spatial transcriptomics, CITE-seq",
        "significance": "Checkpoint inhibitor response biomarker; single-cell resolution "
                        "reveals PD-L1 expression on tumor cells vs. immune cells "
                        "(macrophages, dendritic cells); spatial analysis shows PD-L1 "
                        "distribution relative to CD8+ T cells",
        "workflows": "tme_classification,drug_response",
    },
    "exhaustion-signature": {
        "full_name": "T Cell Exhaustion Gene Signature",
        "assay": "scRNA-seq (PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, TOX, ENTPD1)",
        "significance": "Multi-gene exhaustion program in tumor-infiltrating CD8+ T cells; "
                        "predicts immunotherapy resistance; TOX expression level correlates "
                        "with exhaustion depth; progenitor exhausted (TCF7+) vs. terminally "
                        "exhausted (TCF7-) distinction informs prognosis",
        "workflows": "tme_classification,treatment_monitoring",
    },
    "stemness-score": {
        "full_name": "Cancer Stemness Score (CytoTRACE / scRNA-seq)",
        "assay": "scRNA-seq with CytoTRACE, stemness gene signatures",
        "significance": "Computational measure of differentiation potential from scRNA-seq; "
                        "high stemness correlates with therapy resistance, metastatic "
                        "potential, and poor prognosis; identifies cancer stem cell "
                        "populations for targeted intervention",
        "workflows": "subclonal_architecture,biomarker_discovery",
    },
    "rna-velocity": {
        "full_name": "RNA Velocity (scVelo / Velocyto)",
        "assay": "scRNA-seq (spliced/unspliced ratio)",
        "significance": "Infers future transcriptional state from spliced/unspliced mRNA "
                        "ratios; reveals differentiation directionality, cell fate "
                        "bifurcations, and dynamic state transitions; identifies cells "
                        "transitioning toward drug-resistant states",
        "workflows": "trajectory_inference,treatment_monitoring",
    },
    "cnv-inference": {
        "full_name": "Copy Number Variation Inference (inferCNV / CopyKAT)",
        "assay": "scRNA-seq computational inference",
        "significance": "Distinguishes malignant from non-malignant cells by inferring "
                        "large-scale chromosomal gains/losses from gene expression; "
                        "identifies subclonal CNV events and clonal architecture; "
                        "validated against bulk WGS/WES",
        "workflows": "subclonal_architecture,cell_type_annotation",
    },
    "spatial-autocorrelation": {
        "full_name": "Spatial Autocorrelation (Moran's I / SpatialDE)",
        "assay": "Spatial transcriptomics (Visium, MERFISH, Slide-seq)",
        "significance": "Identifies spatially variable genes beyond random distribution; "
                        "Moran's I > 0 indicates spatial clustering; SpatialDE identifies "
                        "genes with significant spatial expression patterns for niche "
                        "identification and tissue architecture characterization",
        "workflows": "spatial_analysis,biomarker_discovery",
    },
    "clonotype-diversity": {
        "full_name": "TCR/BCR Clonotype Diversity",
        "assay": "scRNA-seq + VDJ sequencing (10x Chromium 5')",
        "significance": "Measures immune repertoire diversity at single-cell resolution; "
                        "clonal expansion of specific TCR clonotypes indicates antigen-driven "
                        "response; shared clonotypes across tissues suggest systemic immune "
                        "response; Shannon entropy and normalized diversity indices",
        "workflows": "treatment_monitoring,biomarker_discovery",
    },
    "ligand-receptor-score": {
        "full_name": "Ligand-Receptor Interaction Score (CellChat / LIANA+)",
        "assay": "scRNA-seq computational inference",
        "significance": "Quantifies intercellular communication probability from "
                        "co-expression of ligand-receptor pairs across cell types; "
                        "identifies signaling pathways driving TME crosstalk, "
                        "paracrine loops, and juxtacrine interactions; prioritizes "
                        "druggable communication axes",
        "workflows": "ligand_receptor,spatial_analysis",
    },
    "cell-cycle-score": {
        "full_name": "Cell Cycle Phase Score (S / G2M / G1)",
        "assay": "scRNA-seq (Seurat/Scanpy cell cycle scoring)",
        "significance": "Assigns cell cycle phase using canonical S and G2M gene sets; "
                        "identifies proliferative tumor subpopulations; confounding factor "
                        "in clustering that may require regression; high proliferation "
                        "index in cancer cells correlates with aggression and chemo-sensitivity",
        "workflows": "cell_type_annotation,subclonal_architecture",
    },
    "gene-module-score": {
        "full_name": "Gene Module / Program Score (NMF / cNMF / Hotspot)",
        "assay": "scRNA-seq non-negative matrix factorization",
        "significance": "Identifies co-expressed gene programs (transcriptional modules) "
                        "active in specific cell populations; discovers tumor-intrinsic "
                        "programs (EMT, hypoxia, stress response, immune evasion) and "
                        "cell-state-specific signatures; consensus NMF (cNMF) provides "
                        "robust program identification across samples",
        "workflows": "biomarker_discovery,subclonal_architecture",
    },
    "adt-protein-level": {
        "full_name": "Antibody-Derived Tag (ADT) Protein Quantification",
        "assay": "CITE-seq / REAP-seq (protein + RNA co-measurement)",
        "significance": "Direct protein quantification at single-cell level using "
                        "oligo-conjugated antibodies; validates RNA-level predictions; "
                        "resolves discordance between mRNA and protein (post-transcriptional "
                        "regulation); critical for surface marker validation in CAR-T and "
                        "immunophenotyping applications",
        "workflows": "cart_target_validation,cell_type_annotation",
    },
    "chromatin-accessibility": {
        "full_name": "Chromatin Accessibility (scATAC-seq peaks)",
        "assay": "scATAC-seq, Multiome (10x), SHARE-seq",
        "significance": "Maps open chromatin regions at single-cell resolution; identifies "
                        "cell-type-specific regulatory elements, transcription factor motifs, "
                        "and gene regulatory networks; Multiome links chromatin state to "
                        "gene expression in the same cell; detects epigenetic reprogramming "
                        "events in cancer and development",
        "workflows": "cell_type_annotation,trajectory_inference",
    },
    "depmap-dependency": {
        "full_name": "DepMap Genetic Dependency Score",
        "assay": "CRISPR screen (Broad DepMap portal)",
        "significance": "Genome-wide genetic dependency scores from cancer cell line CRISPR "
                        "screens; identifies essential genes in specific cancer types; "
                        "cross-referenced with scRNA-seq to nominate cell-type-specific "
                        "therapeutic targets; negative dependency scores indicate essentiality",
        "workflows": "drug_response,biomarker_discovery",
    },
    "immune-exclusion-score": {
        "full_name": "Immune Exclusion Score (TGF-beta + Collagen Density)",
        "assay": "spatial transcriptomics, scRNA-seq with deconvolution",
        "significance": "Composite metric of TGF-beta signaling activity and collagen gene "
                        "expression (COL1A1, COL1A2, COL3A1) at the tumor-stroma interface; "
                        "high scores indicate fibrotic barrier preventing T cell infiltration; "
                        "predicts anti-PD-1 resistance in excluded TME phenotypes",
        "workflows": "tme_classification,spatial_analysis",
    },
    "t-cell-clonality": {
        "full_name": "T Cell Clonality (TCR Repertoire Diversity)",
        "assay": "scTCR-seq, 10x 5' VDJ sequencing",
        "significance": "Shannon entropy-based measure of TCR repertoire diversity; low "
                        "entropy (high clonality) indicates antigen-driven clonal expansion; "
                        "monoclonal expansion in TILs correlates with neoantigen recognition; "
                        "clonality shifts under immunotherapy track treatment response",
        "workflows": "treatment_monitoring,biomarker_discovery",
    },
    "m1-m2-ratio": {
        "full_name": "M1/M2 Macrophage Polarization Ratio",
        "assay": "scRNA-seq (gene signature scoring)",
        "significance": "Ratio of pro-inflammatory M1 (TNF, IL1B, NOS2, CD80) to "
                        "anti-inflammatory M2 (CD163, MRC1, MSR1, IL10) macrophage gene "
                        "signatures; high M1/M2 ratio correlates with immunotherapy response; "
                        "spatial M1/M2 distribution reveals functional TME heterogeneity",
        "workflows": "tme_classification,drug_response",
    },
    "cancer-stemness": {
        "full_name": "Cancer Stemness Signature (ALDH1A1 / CD44 / SOX2)",
        "assay": "scRNA-seq, CITE-seq",
        "significance": "Multi-marker stemness program combining ALDH1A1, CD44, SOX2, NANOG, "
                        "and OCT4 expression; identifies tumor-initiating cell populations; "
                        "high stemness correlates with chemoresistance, metastatic capacity, "
                        "and poor overall survival across cancer types",
        "workflows": "subclonal_architecture,drug_response",
    },
    "metabolic-fitness": {
        "full_name": "Metabolic Fitness Score (OXPHOS vs Glycolysis Ratio)",
        "assay": "scRNA-seq (metabolic pathway gene scoring)",
        "significance": "Ratio of oxidative phosphorylation (OXPHOS) to glycolysis gene "
                        "signatures in immune and tumor cells; metabolically fit T cells "
                        "(high OXPHOS) show superior anti-tumor activity; tumor cells with "
                        "high glycolysis (Warburg effect) create immunosuppressive niche",
        "workflows": "tme_classification,biomarker_discovery",
    },
    "spatial-colocalization": {
        "full_name": "Spatial Colocalization Score (Tumor-Immune Proximity)",
        "assay": "spatial transcriptomics, CODEX, MERFISH",
        "significance": "Quantifies physical proximity between tumor cells and immune cell "
                        "types using Ripley's K, neighborhood enrichment, or nearest-neighbor "
                        "distance metrics; CD8+ T cell colocalization with tumor cells predicts "
                        "immunotherapy response; spatial exclusion patterns identify cold niches",
        "workflows": "spatial_analysis,tme_classification",
    },
    "antigen-presentation": {
        "full_name": "Antigen Presentation Machinery Score (HLA-I + B2M + TAP1)",
        "assay": "scRNA-seq, spatial transcriptomics",
        "significance": "Composite expression of MHC class I components (HLA-A, HLA-B, HLA-C), "
                        "beta-2-microglobulin (B2M), and peptide transporter (TAP1, TAP2); "
                        "downregulation indicates immune evasion mechanism; B2M loss is a "
                        "common resistance mechanism to checkpoint immunotherapy",
        "workflows": "tme_classification,drug_response",
    },
    "proliferation-index": {
        "full_name": "Proliferation Index (MKI67 + TOP2A Composite)",
        "assay": "scRNA-seq, spatial transcriptomics, immunohistochemistry",
        "significance": "Composite score of MKI67 (Ki-67) and TOP2A expression at single-cell "
                        "resolution; identifies actively cycling tumor subpopulations; high "
                        "proliferation index correlates with chemosensitivity but aggressive "
                        "biology; spatial mapping reveals proliferative tumor niches",
        "workflows": "biomarker_discovery,subclonal_architecture",
    },
}


# =====================================================================
# SEARCH PLAN DATACLASS
# =====================================================================

@dataclass
class SearchPlan:
    """Agent's plan for answering a single-cell intelligence question.

    The search plan captures all entities detected in the user's question
    and the strategy the agent will use to retrieve evidence from the
    12 single-cell-specific Milvus collections.
    """
    question: str
    conditions: List[str] = field(default_factory=list)
    cell_types: List[str] = field(default_factory=list)
    biomarkers: List[str] = field(default_factory=list)
    relevant_workflows: List[SCWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, differential, exploratory
    sub_questions: List[str] = field(default_factory=list)
    identified_topics: List[str] = field(default_factory=list)


# =====================================================================
# SINGLE-CELL INTELLIGENCE AGENT
# =====================================================================

class SingleCellAgent:
    """Autonomous Single-Cell Intelligence Agent.

    Wraps the multi-collection SingleCellRAGEngine with planning and reasoning
    capabilities. Designed to answer complex cross-functional questions
    about single-cell genomics, spatial biology, and translational applications.

    Example queries this agent handles:
    - "Annotate cell types in a NSCLC scRNA-seq dataset with 15 clusters"
    - "Classify the tumor microenvironment as hot, cold, excluded, or immunosuppressive"
    - "Predict drug response to pembrolizumab using scRNA-seq TME features"
    - "Detect subclonal architecture in AML from scRNA-seq with inferred CNVs"
    - "Identify spatial niches in a Visium breast cancer dataset"
    - "Infer differentiation trajectory from HSC to mature myeloid cells"
    - "Map ligand-receptor interactions between CAFs and tumor cells"
    - "Validate CD19 as CAR-T target -- check on-tumor vs off-tumor expression"
    - "Compare pre-treatment and post-treatment scRNA-seq for therapy response"
    - "Discover cell-type-specific biomarkers for immunotherapy response prediction"

    Usage:
        agent = SingleCellAgent(rag_engine)
        plan = agent.search_plan("Classify TME in NSCLC scRNA-seq with spatial data")
        response = agent.run("Classify TME in NSCLC scRNA-seq with spatial data")
    """

    def __init__(self, rag_engine):
        """Initialize agent with a configured RAG engine.

        Args:
            rag_engine: SingleCellRAGEngine instance with Milvus collections connected.
        """
        self.rag = rag_engine
        self.knowledge = {
            "conditions": SC_CONDITIONS,
            "cell_types": SC_CELL_TYPES,
            "biomarkers": SC_BIOMARKERS,
        }

    # -- Public API --------------------------------------------------------

    def run(
        self,
        query: str,
        workflow_type: Optional[SCWorkflowType] = None,
        patient_context: Optional[dict] = None,
        **kwargs,
    ) -> SCResponse:
        """Execute the full agent pipeline: plan -> search -> evaluate -> synthesize.

        Args:
            query: Natural language question about single-cell genomics.
            workflow_type: Optional workflow override for collection boosting.
            patient_context: Optional patient/sample data for context.
            **kwargs: Additional query parameters (top_k, collection_filter).

        Returns:
            SCResponse with findings, recommendations, and metadata.
        """
        # Phase 1: Plan
        plan = self.search_plan(query)

        # Phase 2: Determine workflow (allow override)
        workflow = workflow_type or (
            plan.relevant_workflows[0] if plan.relevant_workflows else None
        )

        # Phase 3: Search via RAG engine
        top_k = kwargs.get("top_k", 5)

        response = self.rag.query(
            question=query,
            workflow=workflow,
            top_k=top_k,
            patient_context=patient_context,
        )

        # Phase 4: Evaluate and potentially expand
        if hasattr(response, "results") and response.results is not None:
            quality = self.evaluate_evidence(response.results)
            if quality == "insufficient" and plan.sub_questions:
                for sub_q in plan.sub_questions[:2]:
                    sub_response = self.rag.search(sub_q, top_k=top_k)
                    if sub_response:
                        response.results.extend(sub_response)

        # Phase 5: Check for clinical alerts
        if hasattr(response, "clinical_alerts"):
            response.clinical_alerts = self._check_clinical_alerts(query, plan)

        return response

    def search_plan(self, question: str) -> SearchPlan:
        """Analyze a question and create an optimised search plan.

        Detects single-cell conditions, cell types, and biomarkers in the
        question text. Determines relevant SC workflows, chooses a search
        strategy, and generates sub-questions for comprehensive retrieval
        across collections.

        Args:
            question: The user's natural language question.

        Returns:
            SearchPlan with all detected entities and retrieval strategy.
        """
        plan = SearchPlan(question=question)

        # Step 1: Detect entities
        entities = self._detect_entities(question)
        plan.conditions = entities.get("conditions", [])
        plan.cell_types = entities.get("cell_types", [])
        plan.biomarkers = entities.get("biomarkers", [])

        # Step 2: Determine relevant workflows
        plan.relevant_workflows = [self._detect_workflow(question)]
        # Add entity-derived workflows
        for condition in plan.conditions:
            info = SC_CONDITIONS.get(condition, {})
            for wf in info.get("workflows", []):
                if wf not in plan.relevant_workflows:
                    plan.relevant_workflows.append(wf)

        # Step 3: Choose search strategy
        plan.search_strategy = self._choose_strategy(
            question, plan.conditions, plan.cell_types,
        )

        # Step 4: Generate sub-questions
        plan.sub_questions = self._generate_sub_questions(plan)

        # Step 5: Compile identified topics
        plan.identified_topics = (
            plan.conditions + plan.cell_types + plan.biomarkers
        )

        return plan

    def evaluate_evidence(self, results) -> str:
        """Evaluate the quality and coverage of retrieved evidence.

        Uses collection diversity and hit count to assess whether
        the retrieved evidence is sufficient for a comprehensive answer.

        Args:
            results: List of search results from the RAG engine.

        Returns:
            "sufficient", "partial", or "insufficient".
        """
        if not results:
            return "insufficient"

        total_hits = len(results)
        collections_seen = set()

        for result in results:
            if hasattr(result, "collection"):
                collections_seen.add(result.collection)
            elif isinstance(result, dict):
                collections_seen.add(result.get("collection", "unknown"))

        num_collections = len(collections_seen)

        if num_collections >= 3 and total_hits >= 10:
            return "sufficient"
        elif num_collections >= 2 and total_hits >= 5:
            return "partial"
        else:
            return "insufficient"

    def generate_report(
        self,
        results,
        workflow: Optional[SCWorkflowType] = None,
    ) -> str:
        """Generate a structured single-cell intelligence report.

        Args:
            results: Response object from run() or rag.query().
            workflow: Optional workflow type for section customisation.

        Returns:
            Formatted markdown report string.
        """
        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        question = results.question if hasattr(results, "question") else ""
        plan = self.search_plan(question) if question else SearchPlan(question="")

        report_lines = [
            "# Single-Cell Intelligence Report",
            f"**Query:** {question}",
            f"**Generated:** {timestamp}",
            f"**Workflows:** {', '.join(wf.value for wf in plan.relevant_workflows)}",
            f"**Strategy:** {plan.search_strategy}",
            "",
        ]

        # Detected entities
        if plan.conditions or plan.cell_types or plan.biomarkers:
            report_lines.extend([
                "---",
                "",
                "## Detected Entities",
                "",
            ])
            if plan.conditions:
                report_lines.append(
                    f"- **Conditions / Tissues:** {', '.join(plan.conditions)}"
                )
            if plan.cell_types:
                report_lines.append(
                    f"- **Cell Types:** {', '.join(plan.cell_types)}"
                )
            if plan.biomarkers:
                report_lines.append(
                    f"- **Biomarkers / Assays:** {', '.join(plan.biomarkers)}"
                )
            report_lines.append("")

        # Clinical alerts check
        alerts = self._check_clinical_alerts(question, plan)
        if alerts:
            report_lines.extend([
                "---",
                "",
                "## [CRITICAL] Alerts",
                "",
            ])
            for alert in alerts:
                report_lines.append(f"- **[CRITICAL]** {alert}")
            report_lines.append("")

        # Critical findings from results
        critical_flags = []
        if hasattr(results, "results") and results.results:
            for r in results.results:
                meta = r.metadata if hasattr(r, "metadata") else {}
                if meta.get("urgency") == "critical" or meta.get("safety_alert"):
                    critical_flags.append(r)

        if critical_flags:
            if not alerts:
                report_lines.extend([
                    "---",
                    "",
                    "## [CRITICAL] Safety Alerts",
                    "",
                ])
            for flag in critical_flags:
                text = flag.text if hasattr(flag, "text") else str(flag)
                report_lines.append(
                    f"- **[CRITICAL]** {text[:200]} -- "
                    f"immediate expert review required."
                )
            report_lines.append("")

        # Analysis section
        report_lines.extend([
            "---",
            "",
            "## Analysis",
            "",
        ])

        if hasattr(results, "answer"):
            report_lines.append(results.answer)
        elif hasattr(results, "summary"):
            report_lines.append(results.summary)
        elif isinstance(results, str):
            report_lines.append(results)
        else:
            report_lines.append("No analysis generated.")

        report_lines.append("")

        # Workflow-specific reference sections
        if workflow == SCWorkflowType.CELL_TYPE_ANNOTATION:
            report_lines.extend([
                "---",
                "",
                "## Reference Databases",
                "",
                "- Human Cell Atlas (HCA): https://www.humancellatlas.org/",
                "- Cell Ontology: https://www.ebi.ac.uk/ols4/ontologies/cl",
                "- CellMarker 2.0: http://bio-bigdata.hrbmu.edu.cn/CellMarker/",
                "- PanglaoDB: https://panglaodb.se/",
                "- CellTypist: https://www.celltypist.org/",
                "",
            ])
        elif workflow == SCWorkflowType.TME_CLASSIFICATION:
            report_lines.extend([
                "---",
                "",
                "## Reference Frameworks",
                "",
                "- Galon & Bruni 2019: Immunoscore and immunoprofiling (PMID:30804515)",
                "- Bagaev et al. 2021: TME subtypes from transcriptomics (PMID:33504936)",
                "- Binnewies et al. 2018: Understanding TME in immunotherapy (PMID:29686425)",
                "- Chen & Mellman 2017: Immune exclusion phenotype (PMID:28187290)",
                "",
            ])
        elif workflow == SCWorkflowType.CART_TARGET_VALIDATION:
            report_lines.extend([
                "---",
                "",
                "## Reference Frameworks",
                "",
                "- Human Protein Atlas: tissue expression for on-target/off-tumor risk",
                "- DepMap: genetic dependency for target essentiality assessment",
                "- CARTography: CAR-T target expression databases",
                "- Tabula Sapiens: multi-organ single-cell atlas for normal tissue expression",
                "",
            ])
        elif workflow == SCWorkflowType.TRAJECTORY_INFERENCE:
            report_lines.extend([
                "---",
                "",
                "## Method References",
                "",
                "- Saelens et al. 2019: Trajectory inference benchmarking (PMID:30936559)",
                "- La Manno et al. 2018: RNA velocity (PMID:30089906)",
                "- Bergen et al. 2020: scVelo generalized RNA velocity (PMID:32747759)",
                "- Lange et al. 2022: CellRank for fate mapping (PMID:35027767)",
                "",
            ])
        elif workflow == SCWorkflowType.SPATIAL_ANALYSIS:
            report_lines.extend([
                "---",
                "",
                "## Method References",
                "",
                "- Kleshchevnikov et al. 2022: cell2location (PMID:35027729)",
                "- Biancalani et al. 2021: Tangram spatial mapping (PMID:34711971)",
                "- Svensson et al. 2018: SpatialDE (PMID:29553579)",
                "- Palla et al. 2022: Squidpy spatial analysis (PMID:35027729)",
                "",
            ])

        # Confidence and metadata
        confidence = results.confidence if hasattr(results, "confidence") else 0.0
        report_lines.extend([
            "---",
            "",
            "## Metadata",
            "",
            f"- **Confidence Score:** {confidence:.3f}",
            f"- **Collections Searched:** {results.collections_searched if hasattr(results, 'collections_searched') else 'N/A'}",
            f"- **Search Time:** {results.search_time_ms if hasattr(results, 'search_time_ms') else 'N/A'} ms",
            "",
            "---",
            "",
            "*This report is generated by the Single-Cell Intelligence Agent "
            "within the HCLS AI Factory. Cell type annotations are computational "
            "predictions requiring experimental validation. Drug response predictions "
            "are hypothesis-generating and require functional confirmation. All findings "
            "should be reviewed by expert bioinformaticians and domain scientists.*",
        ])

        return "\n".join(report_lines)

    # -- Clinical Alert Detection ----------------------------------------

    def _check_clinical_alerts(
        self,
        question: str,
        plan: SearchPlan,
    ) -> List[str]:
        """Check for findings requiring urgent alerts.

        Scans the query text and detected entities for patterns that
        indicate safety-critical single-cell findings.

        Args:
            question: Original query text.
            plan: SearchPlan with detected entities.

        Returns:
            List of alert strings for critical findings.
        """
        alerts: List[str] = []
        text_upper = question.upper()

        # CAR-T on-tumor/off-tumor safety concern
        cart_keywords = ["CAR-T", "CAR T", "CHIMERIC ANTIGEN", "CART"]
        off_tumor_keywords = ["OFF-TUMOR", "OFF TUMOR", "NORMAL TISSUE",
                              "ON-TARGET OFF-TUMOR", "SAFETY", "TOXICITY"]
        if any(kw in text_upper for kw in cart_keywords):
            if any(kw in text_upper for kw in off_tumor_keywords):
                alerts.append(
                    "CAR-T target safety evaluation detected. Assess target "
                    "expression across normal tissues using Tabula Sapiens, "
                    "Human Protein Atlas, and GTEx. Flag any expression above "
                    "background in critical organs (heart, lung, CNS, liver). "
                    "Cross-reference with DepMap dependency scores."
                )

        # Immune evasion / MHC loss
        evasion_keywords = ["MHC-I LOSS", "MHC LOSS", "HLA LOSS",
                            "B2M LOSS", "IMMUNE EVASION", "ANTIGEN ESCAPE",
                            "IMMUNE ESCAPE"]
        if any(kw in text_upper for kw in evasion_keywords):
            alerts.append(
                "Immune evasion mechanism detected. MHC-I loss or antigen "
                "escape may render immunotherapy ineffective. Evaluate HLA "
                "allele-specific expression, B2M status, and alternative "
                "antigen presentation pathways. Consider combination "
                "strategies or alternative targets."
            )

        # Drug resistance clone detection
        resistance_keywords = ["DRUG RESISTANCE", "RESISTANT CLONE",
                               "TREATMENT RESISTANCE", "THERAPY RESISTANCE",
                               "RESISTANT SUBCLONE", "RESISTANT POPULATION"]
        if any(kw in text_upper for kw in resistance_keywords):
            alerts.append(
                "Drug-resistant subclone detected or suspected. Evaluate "
                "clonal frequency, resistance mechanism (target mutation, "
                "bypass pathway, lineage switch), and actionable alternatives. "
                "Longitudinal monitoring recommended to track clonal dynamics."
            )

        # Spatial immune exclusion
        exclusion_keywords = ["IMMUNE EXCLUSION", "EXCLUDED PHENOTYPE",
                              "T CELL EXCLUSION", "IMMUNE DESERT",
                              "COLD TUMOR"]
        if any(kw in text_upper for kw in exclusion_keywords):
            alerts.append(
                "Immune-excluded or cold tumor microenvironment detected. "
                "Checkpoint inhibitor monotherapy may have limited efficacy. "
                "Evaluate TGF-beta pathway, WNT/beta-catenin signaling, "
                "and stromal barrier composition. Consider combination "
                "with TME-modifying agents (anti-TGF-beta, oncolytic virus, "
                "radiation)."
            )

        # Malignant transformation
        if any(kw in text_upper for kw in [
            "MALIGNANT TRANSFORMATION", "DYSPLASIA", "PRECANCEROUS",
            "LINEAGE INFIDELITY", "TRANS-DIFFERENTIATION",
        ]):
            alerts.append(
                "Potential malignant transformation or lineage infidelity "
                "detected. Confirm with CNV inference (inferCNV/CopyKAT), "
                "assess clonal architecture, and evaluate tumor suppressor "
                "status. Immediate pathology correlation recommended."
            )

        # Cytokine release syndrome risk
        if any(kw in text_upper for kw in [
            "CYTOKINE STORM", "CRS", "CYTOKINE RELEASE",
            "NEUROTOXICITY", "ICANS",
        ]):
            if any(kw in text_upper for kw in cart_keywords):
                alerts.append(
                    "CAR-T-associated cytokine release syndrome (CRS) or "
                    "immune effector cell-associated neurotoxicity (ICANS) risk. "
                    "Evaluate monocyte/macrophage activation signatures, IL-6 "
                    "pathway, and GM-CSF expression in scRNA-seq data. "
                    "Pre-treatment TME composition may predict CRS severity."
                )

        return alerts

    # -- Workflow Detection -----------------------------------------------

    def _detect_workflow(self, question: str) -> SCWorkflowType:
        """Detect the most relevant workflow from a question.

        Uses keyword-based heuristics to identify which of the 11
        single-cell workflows is most relevant to the query.

        Args:
            question: The user's natural language question.

        Returns:
            Most relevant SCWorkflowType.
        """
        text_upper = question.upper()

        workflow_scores: Dict[SCWorkflowType, float] = {}

        keyword_workflow_map = {
            SCWorkflowType.CELL_TYPE_ANNOTATION: [
                "CELL TYPE", "ANNOTATION", "ANNOTATE", "CELL IDENTITY",
                "CLUSTER IDENTITY", "CELL POPULATION", "CELLTYPE",
                "CELLMARKER", "PANGLAODB", "CELLTYPIST", "AZIMUTH",
                "SINGLER", "CELL ONTOLOGY", "CL:", "SCARCHES",
                "REFERENCE MAPPING", "LABEL TRANSFER", "MARKER GENE",
                "CANONICAL MARKER", "CLUSTER ANNOTATION",
            ],
            SCWorkflowType.TME_CLASSIFICATION: [
                "TUMOR MICROENVIRONMENT", "TME", "IMMUNE INFILTRATE",
                "HOT TUMOR", "COLD TUMOR", "IMMUNE EXCLUDED",
                "IMMUNOSUPPRESSIVE", "IMMUNE INFLAMED",
                "IMMUNE DESERT", "IMMUNOSCORE", "TUMOR IMMUNE",
                "IMMUNE LANDSCAPE", "IMMUNE CONTEXTURE",
                "IMMUNE CELL COMPOSITION", "IMMUNE PHENOTYPE",
                "TREG INFILTRATION", "MDSC", "M2 MACROPHAGE",
                "TUMOR-ASSOCIATED MACROPHAGE", "TAM",
            ],
            SCWorkflowType.DRUG_RESPONSE: [
                "DRUG RESPONSE", "DRUG SENSITIVITY", "DRUG RESISTANCE",
                "THERAPY RESPONSE", "TREATMENT RESPONSE",
                "PHARMACOGENOMICS", "IC50", "DEPMAP",
                "CHEMOSENSITIVITY", "IMMUNOTHERAPY RESPONSE",
                "CHECKPOINT INHIBITOR", "PEMBROLIZUMAB RESPONSE",
                "RESISTANCE MECHANISM", "THERAPEUTIC TARGET",
                "ACTIONABLE", "DRUGGABLE",
            ],
            SCWorkflowType.SUBCLONAL_ARCHITECTURE: [
                "SUBCLONE", "SUBCLONAL", "CLONAL ARCHITECTURE",
                "CLONAL EVOLUTION", "CLONAL HIERARCHY",
                "CLONAL DYNAMICS", "TUMOR HETEROGENEITY",
                "INTRATUMOR HETEROGENEITY", "ITH",
                "INFERCNV", "COPYBAT", "CNV INFERENCE",
                "FOUNDER CLONE", "BRANCHING EVOLUTION",
                "CANCER STEM CELL", "TUMOR INITIATING",
            ],
            SCWorkflowType.SPATIAL_ANALYSIS: [
                "SPATIAL", "VISIUM", "MERFISH", "SLIDE-SEQ",
                "SPATIAL TRANSCRIPTOMICS", "NICHE", "TISSUE ARCHITECTURE",
                "SPATIAL PATTERN", "CELL2LOCATION", "TANGRAM",
                "SQUIDPY", "SPATIALDE", "SPARK",
                "SPATIAL VARIABLE", "SPATIAL DOMAIN",
                "SPATIAL NICHE", "CELLCHARTER", "COLOCALIZATION",
                "NEIGHBORHOOD", "TISSUE REGION",
            ],
            SCWorkflowType.TRAJECTORY_INFERENCE: [
                "TRAJECTORY", "PSEUDOTIME", "LINEAGE",
                "DIFFERENTIATION", "CELL FATE", "RNA VELOCITY",
                "SCVELO", "VELOCYTO", "MONOCLE", "PAGA",
                "CELLRANK", "DEVELOPMENTAL", "BIFURCATION",
                "BRANCHING", "CYTOTRACE", "STEMNESS",
                "PROGENITOR", "MATURATION", "LINEAGE TRACING",
            ],
            SCWorkflowType.LIGAND_RECEPTOR: [
                "LIGAND-RECEPTOR", "LIGAND RECEPTOR", "CELL COMMUNICATION",
                "CELLCHAT", "LIANA", "NICHENET", "COMMOT",
                "INTERCELLULAR", "SIGNALING", "PARACRINE",
                "JUXTACRINE", "RECEPTOR-LIGAND", "CELL-CELL INTERACTION",
                "CELL CROSSTALK", "SIGNALING PATHWAY",
                "CYTOKINE SIGNALING", "CHEMOKINE",
            ],
            SCWorkflowType.BIOMARKER_DISCOVERY: [
                "BIOMARKER", "DIFFERENTIAL EXPRESSION", "DEG",
                "MARKER GENE", "SIGNATURE", "GENE SIGNATURE",
                "PROGNOSTIC", "PREDICTIVE MARKER",
                "CELL-TYPE-SPECIFIC MARKER", "SURFACE MARKER",
                "DIAGNOSTIC MARKER", "GENE MODULE",
                "TRANSCRIPTIONAL PROGRAM", "NMF", "CNMF",
                "HOTSPOT", "GENE SET",
            ],
            SCWorkflowType.CART_TARGET_VALIDATION: [
                "CAR-T", "CAR T", "CHIMERIC ANTIGEN", "CART",
                "CD19 TARGET", "BCMA TARGET", "CD22 TARGET",
                "ON-TUMOR", "OFF-TUMOR", "TARGET VALIDATION",
                "ANTIGEN EXPRESSION", "ANTIGEN ESCAPE",
                "TARGET EXPRESSION", "BISPECIFIC",
                "ANTIBODY-DRUG CONJUGATE", "ADC",
                "GPRC5D", "FCRH5",
            ],
            SCWorkflowType.TREATMENT_MONITORING: [
                "TREATMENT MONITORING", "LONGITUDINAL",
                "PRE-TREATMENT", "POST-TREATMENT", "ON-TREATMENT",
                "RESPONSE MONITORING", "CLONAL DYNAMICS",
                "MRD", "MINIMAL RESIDUAL DISEASE",
                "RELAPSE", "REMISSION", "TREATMENT RESPONSE",
                "BEFORE AND AFTER", "TIMEPOINT",
                "SERIAL SAMPLING", "THERAPY MONITORING",
            ],
        }

        for wf, keywords in keyword_workflow_map.items():
            for kw in keywords:
                if kw in text_upper:
                    workflow_scores[wf] = workflow_scores.get(wf, 0) + 1.0

        if not workflow_scores:
            return SCWorkflowType.GENERAL

        sorted_workflows = sorted(
            workflow_scores.items(), key=lambda x: x[1], reverse=True,
        )

        return sorted_workflows[0][0]

    # -- Entity Detection -------------------------------------------------

    def _detect_entities(self, question: str) -> Dict[str, List[str]]:
        """Detect single-cell entities in the question text.

        Scans for conditions, cell types, and biomarkers using the knowledge
        dictionaries. Performs case-insensitive matching against canonical
        names and aliases.

        Args:
            question: The user's natural language question.

        Returns:
            Dict with keys 'conditions', 'cell_types', 'biomarkers' mapping
            to lists of detected entity names.
        """
        import re

        entities: Dict[str, List[str]] = {
            "conditions": [],
            "cell_types": [],
            "biomarkers": [],
        }

        text_lower = question.lower()

        # Detect conditions
        for condition, info in SC_CONDITIONS.items():
            if condition in text_lower:
                if condition not in entities["conditions"]:
                    entities["conditions"].append(condition)
                continue
            aliases = info.get("aliases", [])
            for alias in aliases:
                if len(alias) <= 3:
                    pattern = r'\b' + re.escape(alias) + r'\b'
                    if re.search(pattern, text_lower):
                        if condition not in entities["conditions"]:
                            entities["conditions"].append(condition)
                        break
                else:
                    if alias.lower() in text_lower:
                        if condition not in entities["conditions"]:
                            entities["conditions"].append(condition)
                        break

        # Detect cell types
        for cell_type, info in SC_CELL_TYPES.items():
            canonical = cell_type.replace("_", " ")
            if canonical in text_lower:
                if cell_type not in entities["cell_types"]:
                    entities["cell_types"].append(cell_type)
                continue
            full_name = info.get("full_name", "")
            if full_name and full_name.lower() in text_lower:
                if cell_type not in entities["cell_types"]:
                    entities["cell_types"].append(cell_type)
                continue
            # Check CL ID
            cl_id = info.get("cell_ontology_id", "")
            if cl_id and cl_id in question:
                if cell_type not in entities["cell_types"]:
                    entities["cell_types"].append(cell_type)
                continue
            # Check markers (if at least 2 canonical markers mentioned)
            markers = info.get("canonical_markers", [])
            markers_found = sum(1 for m in markers if m.upper() in question.upper())
            if markers_found >= 2:
                if cell_type not in entities["cell_types"]:
                    entities["cell_types"].append(cell_type)

        # Detect biomarkers
        for biomarker, info in SC_BIOMARKERS.items():
            canonical = biomarker.replace("-", " ").replace("_", " ")
            if canonical in text_lower:
                if biomarker not in entities["biomarkers"]:
                    entities["biomarkers"].append(biomarker)
                continue
            full_name = info.get("full_name", "")
            if full_name and full_name.lower() in text_lower:
                if biomarker not in entities["biomarkers"]:
                    entities["biomarkers"].append(biomarker)

        return entities

    # -- Search Strategy ---------------------------------------------------

    def _choose_strategy(
        self,
        text: str,
        conditions: List[str],
        cell_types: List[str],
    ) -> str:
        """Choose search strategy: broad, targeted, differential, or exploratory.

        Args:
            text: Original query text.
            conditions: Detected conditions.
            cell_types: Detected cell types.

        Returns:
            Strategy name string.
        """
        text_upper = text.upper()

        # Exploratory queries (discovery-oriented)
        exploratory_keywords = [
            "DISCOVER", "NOVEL", "UNKNOWN", "UNCHARACTERIZED",
            "DE NOVO", "EXPLORE", "WHAT ARE", "IDENTIFY NEW",
            "FIND NEW", "CHARACTERIZE",
        ]
        if any(kw in text_upper for kw in exploratory_keywords):
            return "exploratory"

        # Differential / comparison queries
        if ("DIFFERENTIAL" in text_upper or "COMPARE" in text_upper
                or "DISTINGUISH" in text_upper or " VS " in text_upper
                or "VERSUS" in text_upper or "DIFFERENCE BETWEEN" in text_upper
                or "PRE-TREATMENT" in text_upper or "POST-TREATMENT" in text_upper):
            return "differential"

        # Targeted: specific condition + cell type or single focused entity
        if (len(conditions) == 1 and len(cell_types) <= 2) or (
            len(conditions) <= 1 and len(cell_types) == 1
        ):
            if conditions or cell_types:
                return "targeted"

        return "broad"

    # -- Sub-Question Generation -------------------------------------------

    def _generate_sub_questions(self, plan: SearchPlan) -> List[str]:
        """Generate sub-questions for comprehensive retrieval.

        Decomposes the main question into focused sub-queries based on
        the detected entities and workflow type. Enables multi-hop
        retrieval across different aspects of the single-cell question.

        Args:
            plan: SearchPlan with detected entities and workflows.

        Returns:
            List of sub-question strings (typically 2-4 questions).
        """
        sub_questions: List[str] = []

        condition_label = plan.conditions[0] if plan.conditions else "the tissue/condition"
        cell_type_label = plan.cell_types[0] if plan.cell_types else "the cell population"
        plan.biomarkers[0] if plan.biomarkers else "the biomarker"

        primary_wf = (
            plan.relevant_workflows[0] if plan.relevant_workflows
            else SCWorkflowType.GENERAL
        )

        # -- Pattern 1: Cell Type Annotation ----------------------------
        if primary_wf == SCWorkflowType.CELL_TYPE_ANNOTATION:
            sub_questions = [
                f"What are the canonical marker genes for cell types in {condition_label}?",
                f"What reference atlases are available for {condition_label} cell type annotation?",
                f"What computational methods are recommended for annotating {cell_type_label}?",
                f"What Cell Ontology IDs and marker databases describe {cell_type_label}?",
            ]

        # -- Pattern 2: TME Classification ------------------------------
        elif primary_wf == SCWorkflowType.TME_CLASSIFICATION:
            sub_questions = [
                f"What immune cell composition defines the TME subtype in {condition_label}?",
                "What gene signatures distinguish hot, cold, excluded, and immunosuppressive TME?",
                f"What spatial patterns characterize immune infiltration in {condition_label}?",
                f"What are the immunotherapy response implications of TME classification in {condition_label}?",
            ]

        # -- Pattern 3: Drug Response -----------------------------------
        elif primary_wf == SCWorkflowType.DRUG_RESPONSE:
            sub_questions = [
                f"What cell-type-specific drug response signatures exist for {condition_label}?",
                f"What resistance mechanisms are observed at single-cell resolution in {condition_label}?",
                f"What DepMap dependency scores are relevant for druggable targets in {condition_label}?",
                f"What scRNA-seq biomarkers predict therapy response in {condition_label}?",
            ]

        # -- Pattern 4: Subclonal Architecture --------------------------
        elif primary_wf == SCWorkflowType.SUBCLONAL_ARCHITECTURE:
            sub_questions = [
                f"What methods detect subclonal architecture from scRNA-seq in {condition_label}?",
                f"What CNV inference approaches (inferCNV, CopyKAT) are validated for {condition_label}?",
                f"What clonal evolution patterns are observed in {condition_label}?",
                f"How does subclonal heterogeneity affect therapy outcomes in {condition_label}?",
            ]

        # -- Pattern 5: Spatial Analysis --------------------------------
        elif primary_wf == SCWorkflowType.SPATIAL_ANALYSIS:
            sub_questions = [
                f"What spatial transcriptomics technologies are suitable for {condition_label}?",
                f"What spatial niches and domains are characteristic of {condition_label}?",
                f"What cell-cell proximity patterns are significant in {condition_label}?",
                f"What computational methods identify spatially variable genes in {condition_label}?",
            ]

        # -- Pattern 6: Trajectory Inference ----------------------------
        elif primary_wf == SCWorkflowType.TRAJECTORY_INFERENCE:
            sub_questions = [
                f"What differentiation trajectory connects {cell_type_label} to mature cell types?",
                f"What RNA velocity analysis reveals about cell state transitions in {condition_label}?",
                f"What are the key transcription factors driving lineage commitment in {cell_type_label}?",
                f"What pseudotime methods are benchmarked for {condition_label}?",
            ]

        # -- Pattern 7: Ligand-Receptor ---------------------------------
        elif primary_wf == SCWorkflowType.LIGAND_RECEPTOR:
            sub_questions = [
                f"What ligand-receptor interactions are active between cell types in {condition_label}?",
                f"What signaling pathways mediate cell-cell communication in {condition_label}?",
                f"What methods (CellChat, LIANA+, NicheNet) are recommended for {condition_label}?",
                f"What druggable ligand-receptor axes exist in {condition_label}?",
            ]

        # -- Pattern 8: Biomarker Discovery -----------------------------
        elif primary_wf == SCWorkflowType.BIOMARKER_DISCOVERY:
            sub_questions = [
                f"What cell-type-specific biomarkers distinguish {condition_label} from normal?",
                f"What gene signatures or modules are associated with {condition_label} prognosis?",
                f"What surface markers enable therapeutic targeting in {condition_label}?",
                f"What transcriptional programs are unique to {cell_type_label} in {condition_label}?",
            ]

        # -- Pattern 9: CAR-T Target Validation -------------------------
        elif primary_wf == SCWorkflowType.CART_TARGET_VALIDATION:
            sub_questions = [
                f"What is the on-tumor expression of the target in {condition_label} at single-cell level?",
                "What is the off-tumor expression of the target across normal tissues?",
                f"What antigen escape mechanisms are observed in {condition_label}?",
                "What DepMap dependency scores indicate target essentiality for tumor survival?",
            ]

        # -- Pattern 10: Treatment Monitoring ---------------------------
        elif primary_wf == SCWorkflowType.TREATMENT_MONITORING:
            sub_questions = [
                f"What cell population changes occur during treatment in {condition_label}?",
                f"What clonal dynamics are observed between pre- and post-treatment in {condition_label}?",
                f"What immune cell state changes indicate therapy response in {condition_label}?",
                f"What scRNA-seq signatures predict relapse or MRD in {condition_label}?",
            ]

        # -- Default ---------------------------------------------------
        else:
            sub_questions = [
                f"What cell types are present in {condition_label} at single-cell resolution?",
                f"What computational methods are recommended for analyzing {condition_label}?",
                f"What is the current single-cell evidence base for {condition_label}?",
            ]

        return sub_questions

    # -- Build Search Strategy Description ---------------------------------

    def _build_search_strategy(
        self,
        entities: Dict[str, List[str]],
        workflow: SCWorkflowType,
    ) -> str:
        """Build a descriptive search strategy based on entities and workflow.

        Args:
            entities: Detected entities dict from _detect_entities.
            workflow: Determined workflow type.

        Returns:
            Strategy description string for logging/debugging.
        """
        parts = [f"Workflow: {workflow.value}"]

        if entities.get("conditions"):
            parts.append(f"Conditions: {', '.join(entities['conditions'])}")
        if entities.get("cell_types"):
            parts.append(f"Cell Types: {', '.join(entities['cell_types'])}")
        if entities.get("biomarkers"):
            parts.append(f"Biomarkers: {', '.join(entities['biomarkers'])}")

        # Determine collection priorities
        boosts = WORKFLOW_COLLECTION_BOOST.get(workflow, {})
        top_collections = sorted(
            boosts.items(), key=lambda x: x[1], reverse=True,
        )[:5]
        if top_collections:
            parts.append(
                "Priority collections: "
                + ", ".join(f"{c}({w:.1f}x)" for c, w in top_collections)
            )

        return " | ".join(parts)
