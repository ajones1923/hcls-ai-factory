"""Pydantic data models for the Single-Cell Intelligence Agent.

Comprehensive enums and models for a single-cell transcriptomics RAG-based
analysis system covering cell type annotation, tumor microenvironment profiling,
drug response prediction, subclonal architecture analysis, spatial niche
identification, trajectory analysis, ligand-receptor interactions,
biomarker discovery, CAR-T target validation, and treatment monitoring.

Follows the same dataclass/Pydantic pattern as:
  - rare_disease_diagnostic_agent/src/models.py
  - clinical_trial_intelligence_agent/src/models.py

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


# ===================================================================
# ENUMS
# ===================================================================


class SCWorkflowType(str, Enum):
    """Types of single-cell analysis query workflows."""
    CELL_TYPE_ANNOTATION = "cell_type_annotation"
    TME_PROFILING = "tme_profiling"
    DRUG_RESPONSE = "drug_response"
    SUBCLONAL_ARCHITECTURE = "subclonal_architecture"
    SPATIAL_NICHE = "spatial_niche"
    TRAJECTORY_ANALYSIS = "trajectory_analysis"
    LIGAND_RECEPTOR = "ligand_receptor"
    BIOMARKER_DISCOVERY = "biomarker_discovery"
    CART_TARGET_VALIDATION = "cart_target_validation"
    TREATMENT_MONITORING = "treatment_monitoring"
    GENERAL = "general"


class TMEClass(str, Enum):
    """Tumor microenvironment immunophenotype classification."""
    HOT_INFLAMED = "hot_inflamed"
    COLD_DESERT = "cold_desert"
    EXCLUDED = "excluded"
    IMMUNOSUPPRESSIVE = "immunosuppressive"


class CellTypeConfidence(str, Enum):
    """Confidence level for cell type annotation."""
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class SpatialPlatform(str, Enum):
    """Spatial transcriptomics technology platforms."""
    VISIUM = "visium"
    MERFISH = "merfish"
    XENIUM = "xenium"
    CODEX = "codex"


class ResistanceRisk(str, Enum):
    """Predicted risk of therapeutic resistance."""
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class SeverityLevel(str, Enum):
    """Clinical finding severity classification."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INFORMATIONAL = "informational"


class EvidenceLevel(str, Enum):
    """Strength of evidence supporting an analytical finding."""
    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"
    CONFLICTING = "conflicting"
    UNCERTAIN = "uncertain"


class TrajectoryType(str, Enum):
    """Types of cellular trajectory or differentiation path."""
    DIFFERENTIATION = "differentiation"
    ACTIVATION = "activation"
    EXHAUSTION = "exhaustion"
    EMT = "emt"
    STEMNESS = "stemness"
    CELL_CYCLE = "cell_cycle"


class AssayType(str, Enum):
    """Single-cell assay modalities."""
    SCRNASEQ = "scRNA-seq"
    SNRNASEQ = "snRNA-seq"
    SCATACSEQ = "scATAC-seq"
    CITE_SEQ = "CITE-seq"
    MULTIOME = "multiome"
    SPATIAL = "spatial"


class NormalizationMethod(str, Enum):
    """Normalization approaches for single-cell count data."""
    LOG_NORMALIZE = "log_normalize"
    SCRAN = "scran"
    SCT = "sctransform"
    PEARSON_RESIDUALS = "pearson_residuals"


class ClusteringMethod(str, Enum):
    """Cell clustering algorithm options."""
    LEIDEN = "leiden"
    LOUVAIN = "louvain"
    KMEANS = "kmeans"
    PHENOGRAPH = "phenograph"


# ===================================================================
# PYDANTIC MODELS - QUERY
# ===================================================================


class SCQuery(BaseModel):
    """Input query for single-cell analysis.

    Captures the user question, optional dataset context, tissue/disease
    information, and workflow routing parameters.
    """
    query: str = Field(
        ...,
        max_length=10000,
        description="Natural language question about single-cell data or analysis",
    )
    patient_id: Optional[str] = Field(
        default=None,
        description="Patient or sample identifier",
    )
    dataset_id: Optional[str] = Field(
        default=None,
        description="Reference dataset identifier (e.g., CellxGene dataset ID)",
    )
    tissue_type: Optional[str] = Field(
        default=None,
        description="Tissue of origin (e.g., 'PBMC', 'lung', 'tumor biopsy')",
    )
    disease_context: Optional[str] = Field(
        default=None,
        description="Disease context (e.g., 'NSCLC', 'AML', 'Crohn disease')",
    )
    assay_type: Optional[AssayType] = Field(
        default=None,
        description="Single-cell assay modality used",
    )
    workflow_type: Optional[SCWorkflowType] = Field(
        default=None,
        description="Specific analysis workflow; auto-detected if omitted",
    )
    genes_of_interest: List[str] = Field(
        default_factory=list,
        description="List of gene symbols relevant to the query",
    )
    cell_types_of_interest: List[str] = Field(
        default_factory=list,
        description="Cell types to focus on (e.g., ['CD8+ T cells', 'Macrophages'])",
    )
    spatial_platform: Optional[SpatialPlatform] = Field(
        default=None,
        description="Spatial platform if spatial data is involved",
    )
    top_k: int = Field(
        default=10, ge=1, le=100,
        description="Number of results to return per collection",
    )
    include_methods: bool = Field(
        default=True,
        description="Whether to include methodological references in the response",
    )


# ===================================================================
# PYDANTIC MODELS - SEARCH RESULT
# ===================================================================


class SCSearchResult(BaseModel):
    """A single search result from a Milvus collection."""
    collection: str = Field(
        ...,
        description="Name of the source Milvus collection",
    )
    record_id: int = Field(
        ...,
        description="Primary key of the matched record",
    )
    score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Cosine similarity score (0-1)",
    )
    content: str = Field(
        ...,
        description="Textual content of the matched record",
    )
    metadata: Dict = Field(
        default_factory=dict,
        description="Additional metadata fields from the collection",
    )


# ===================================================================
# PYDANTIC MODELS - CELL TYPE ANNOTATION
# ===================================================================


class CellTypeAnnotation(BaseModel):
    """Result of cell type annotation for a cluster or cell population."""
    cluster_id: str = Field(
        ...,
        description="Cluster identifier from the analysis",
    )
    cell_type: str = Field(
        ...,
        description="Predicted cell type label (e.g., 'CD8+ effector T cell')",
    )
    cell_ontology_id: Optional[str] = Field(
        default=None,
        description="Cell Ontology (CL) identifier",
    )
    confidence: CellTypeConfidence = Field(
        ...,
        description="Confidence level of the annotation",
    )
    confidence_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Numeric confidence score",
    )
    marker_genes: List[str] = Field(
        default_factory=list,
        description="Key marker genes supporting this annotation",
    )
    marker_evidence: Dict[str, float] = Field(
        default_factory=dict,
        description="Gene-to-expression-score mapping for top markers",
    )
    cell_count: int = Field(
        default=0, ge=0,
        description="Number of cells assigned to this type",
    )
    fraction: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Fraction of total cells in the dataset",
    )
    alternative_labels: List[str] = Field(
        default_factory=list,
        description="Alternative cell type labels considered",
    )
    reference_dataset: Optional[str] = Field(
        default=None,
        description="Reference atlas or dataset used for annotation",
    )


# ===================================================================
# PYDANTIC MODELS - TME PROFILE
# ===================================================================


class TMEProfile(BaseModel):
    """Tumor microenvironment profile from single-cell analysis."""
    tme_class: TMEClass = Field(
        ...,
        description="Overall TME immunophenotype classification",
    )
    immune_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Aggregate immune infiltration score",
    )
    stromal_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Aggregate stromal component score",
    )
    cell_type_fractions: Dict[str, float] = Field(
        default_factory=dict,
        description="Fraction of each cell type in the TME",
    )
    exhaustion_signature: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="T cell exhaustion signature score",
    )
    checkpoint_expression: Dict[str, float] = Field(
        default_factory=dict,
        description="Expression levels of immune checkpoint genes (e.g., PD-L1, CTLA-4)",
    )
    cytokine_milieu: Dict[str, float] = Field(
        default_factory=dict,
        description="Key cytokine expression levels in the TME",
    )
    predicted_immunotherapy_response: Optional[str] = Field(
        default=None,
        description="Predicted immunotherapy response (responder/non-responder/uncertain)",
    )
    spatial_pattern: Optional[str] = Field(
        default=None,
        description="Spatial distribution pattern of immune cells",
    )
    evidence_level: EvidenceLevel = Field(
        default=EvidenceLevel.UNCERTAIN,
        description="Strength of evidence for the TME classification",
    )


# ===================================================================
# PYDANTIC MODELS - DRUG RESPONSE PREDICTION
# ===================================================================


class DrugResponsePrediction(BaseModel):
    """Predicted drug response based on single-cell transcriptomic signatures."""
    drug_name: str = Field(
        ...,
        description="Name of the drug or compound",
    )
    drug_class: Optional[str] = Field(
        default=None,
        description="Drug class (e.g., 'checkpoint inhibitor', 'TKI')",
    )
    predicted_sensitivity: float = Field(
        ..., ge=0.0, le=1.0,
        description="Predicted sensitivity score (higher = more sensitive)",
    )
    resistance_risk: ResistanceRisk = Field(
        ...,
        description="Predicted risk of resistance",
    )
    resistance_mechanisms: List[str] = Field(
        default_factory=list,
        description="Identified potential resistance mechanisms",
    )
    resistant_subpopulation: Optional[str] = Field(
        default=None,
        description="Cell subpopulation most likely to confer resistance",
    )
    resistant_fraction: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Estimated fraction of resistant cells",
    )
    synergy_candidates: List[str] = Field(
        default_factory=list,
        description="Drugs that may synergize based on subpopulation vulnerabilities",
    )
    evidence_level: EvidenceLevel = Field(
        default=EvidenceLevel.UNCERTAIN,
        description="Strength of evidence for the prediction",
    )
    source_studies: List[str] = Field(
        default_factory=list,
        description="Key studies supporting this prediction",
    )


# ===================================================================
# PYDANTIC MODELS - SUBCLONAL ARCHITECTURE
# ===================================================================


class SubclonalResult(BaseModel):
    """Subclonal architecture analysis from single-cell data."""
    clone_id: str = Field(
        ...,
        description="Clone identifier",
    )
    clone_fraction: float = Field(
        ..., ge=0.0, le=1.0,
        description="Fraction of cells belonging to this clone",
    )
    cell_count: int = Field(
        default=0, ge=0,
        description="Number of cells in this clone",
    )
    driver_mutations: List[str] = Field(
        default_factory=list,
        description="Driver mutations characterizing this clone",
    )
    cnv_profile: Dict[str, str] = Field(
        default_factory=dict,
        description="Copy number variation profile (chr region -> gain/loss/neutral)",
    )
    transcriptomic_signature: List[str] = Field(
        default_factory=list,
        description="Top differentially expressed genes for this clone",
    )
    phylogenetic_position: Optional[str] = Field(
        default=None,
        description="Position in the clonal phylogeny (e.g., 'root', 'branch_A')",
    )
    fitness_score: float = Field(
        default=0.0, ge=0.0,
        description="Estimated clonal fitness score",
    )
    therapy_implications: List[str] = Field(
        default_factory=list,
        description="Therapeutic implications of this clone's profile",
    )


# ===================================================================
# PYDANTIC MODELS - SPATIAL NICHE
# ===================================================================


class SpatialNiche(BaseModel):
    """Spatial niche identified from spatial transcriptomics data."""
    niche_id: str = Field(
        ...,
        description="Spatial niche identifier",
    )
    niche_label: str = Field(
        ...,
        description="Descriptive label (e.g., 'Tumor-immune interface')",
    )
    dominant_cell_types: List[str] = Field(
        default_factory=list,
        description="Cell types enriched in this spatial niche",
    )
    cell_type_proportions: Dict[str, float] = Field(
        default_factory=dict,
        description="Proportion of each cell type within the niche",
    )
    spatial_platform: Optional[SpatialPlatform] = Field(
        default=None,
        description="Spatial platform used for niche identification",
    )
    signature_genes: List[str] = Field(
        default_factory=list,
        description="Genes characterizing this spatial niche",
    )
    interaction_partners: Dict[str, float] = Field(
        default_factory=dict,
        description="Ligand-receptor interactions enriched in this niche",
    )
    area_fraction: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Fraction of tissue area occupied by this niche",
    )
    clinical_relevance: Optional[str] = Field(
        default=None,
        description="Clinical significance of this spatial niche",
    )


# ===================================================================
# PYDANTIC MODELS - TRAJECTORY RESULT
# ===================================================================


class TrajectoryResult(BaseModel):
    """Cellular trajectory or pseudotime analysis result."""
    trajectory_id: str = Field(
        ...,
        description="Trajectory identifier",
    )
    trajectory_type: TrajectoryType = Field(
        ...,
        description="Type of cellular trajectory",
    )
    start_cell_type: str = Field(
        ...,
        description="Starting cell state of the trajectory",
    )
    end_cell_type: str = Field(
        ...,
        description="Terminal cell state of the trajectory",
    )
    intermediate_states: List[str] = Field(
        default_factory=list,
        description="Intermediate cell states along the trajectory",
    )
    driver_genes: List[str] = Field(
        default_factory=list,
        description="Genes driving the transition along this trajectory",
    )
    branching_points: List[str] = Field(
        default_factory=list,
        description="Points where the trajectory bifurcates",
    )
    pseudotime_range: List[float] = Field(
        default_factory=list,
        description="Pseudotime range [start, end] for this trajectory",
    )
    cell_count: int = Field(
        default=0, ge=0,
        description="Number of cells along this trajectory",
    )
    clinical_relevance: Optional[str] = Field(
        default=None,
        description="Clinical significance of this trajectory",
    )


# ===================================================================
# PYDANTIC MODELS - LIGAND-RECEPTOR INTERACTION
# ===================================================================


class LigandReceptorInteraction(BaseModel):
    """Cell-cell communication via ligand-receptor interaction."""
    ligand_gene: str = Field(
        ...,
        description="Ligand gene symbol",
    )
    receptor_gene: str = Field(
        ...,
        description="Receptor gene symbol",
    )
    source_cell_type: str = Field(
        ...,
        description="Cell type expressing the ligand",
    )
    target_cell_type: str = Field(
        ...,
        description="Cell type expressing the receptor",
    )
    interaction_score: float = Field(
        ..., ge=0.0,
        description="Interaction strength score (higher = stronger)",
    )
    p_value: float = Field(
        default=1.0, ge=0.0, le=1.0,
        description="Statistical significance of the interaction",
    )
    pathway: Optional[str] = Field(
        default=None,
        description="Signaling pathway involved",
    )
    is_novel: bool = Field(
        default=False,
        description="Whether this is a novel (not previously reported) interaction",
    )
    clinical_relevance: Optional[str] = Field(
        default=None,
        description="Therapeutic or diagnostic relevance",
    )
    method: Optional[str] = Field(
        default=None,
        description="Method used to infer interaction (e.g., CellChat, CellPhoneDB)",
    )


# ===================================================================
# PYDANTIC MODELS - BIOMARKER CANDIDATE
# ===================================================================


class BiomarkerCandidate(BaseModel):
    """Biomarker candidate identified from single-cell analysis."""
    gene: str = Field(
        ...,
        description="Gene symbol of the biomarker candidate",
    )
    biomarker_type: str = Field(
        ...,
        description="Type of biomarker (diagnostic, prognostic, predictive, pharmacodynamic)",
    )
    cell_type_specific: str = Field(
        ...,
        description="Cell type in which the biomarker is expressed",
    )
    fold_change: float = Field(
        ...,
        description="Log2 fold change in expression",
    )
    p_value_adjusted: float = Field(
        default=1.0, ge=0.0, le=1.0,
        description="Adjusted p-value (FDR)",
    )
    specificity_score: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Cell-type specificity score (AUROC or similar)",
    )
    surface_protein: bool = Field(
        default=False,
        description="Whether this gene encodes a surface protein (targetable)",
    )
    existing_assay: bool = Field(
        default=False,
        description="Whether a clinical-grade assay exists for this biomarker",
    )
    evidence_level: EvidenceLevel = Field(
        default=EvidenceLevel.UNCERTAIN,
        description="Strength of evidence supporting this candidate",
    )
    supporting_datasets: List[str] = Field(
        default_factory=list,
        description="Datasets in which this biomarker was validated",
    )


# ===================================================================
# PYDANTIC MODELS - CAR-T TARGET VALIDATION
# ===================================================================


class CARTTargetValidation(BaseModel):
    """CAR-T therapy target validation from single-cell expression data."""
    target_gene: str = Field(
        ...,
        description="Target antigen gene symbol",
    )
    tumor_expression_fraction: float = Field(
        ..., ge=0.0, le=1.0,
        description="Fraction of tumor cells expressing the target",
    )
    tumor_expression_level: float = Field(
        ..., ge=0.0,
        description="Mean expression level in tumor cells",
    )
    normal_tissue_expression: Dict[str, float] = Field(
        default_factory=dict,
        description="Expression levels in normal tissue cell types (safety profile)",
    )
    on_target_off_tumor_risk: str = Field(
        default="unknown",
        description="Risk of on-target off-tumor toxicity (high/medium/low/unknown)",
    )
    co_expression_partners: List[str] = Field(
        default_factory=list,
        description="Genes co-expressed with the target for dual-targeting strategies",
    )
    heterogeneity_score: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Expression heterogeneity across tumor subpopulations",
    )
    escape_variant_risk: ResistanceRisk = Field(
        default=ResistanceRisk.MEDIUM,
        description="Risk of antigen-loss escape variants",
    )
    evidence_level: EvidenceLevel = Field(
        default=EvidenceLevel.UNCERTAIN,
        description="Strength of evidence for this target",
    )
    clinical_trials: List[str] = Field(
        default_factory=list,
        description="Active clinical trials targeting this antigen",
    )


# ===================================================================
# PYDANTIC MODELS - TREATMENT RESPONSE
# ===================================================================


class TreatmentResponse(BaseModel):
    """Treatment monitoring result from longitudinal single-cell data."""
    timepoint: str = Field(
        ...,
        description="Treatment timepoint (e.g., 'baseline', 'week_4', 'progression')",
    )
    treatment: str = Field(
        ...,
        description="Treatment regimen being monitored",
    )
    response_category: str = Field(
        default="unknown",
        description="Response category (complete, partial, stable, progression, unknown)",
    )
    compositional_shifts: Dict[str, float] = Field(
        default_factory=dict,
        description="Changes in cell type proportions vs baseline",
    )
    emerging_clones: List[str] = Field(
        default_factory=list,
        description="New clonal populations emerging under treatment",
    )
    resistance_signatures: List[str] = Field(
        default_factory=list,
        description="Gene signatures associated with emerging resistance",
    )
    immune_dynamics: Dict[str, float] = Field(
        default_factory=dict,
        description="Changes in immune cell activation/exhaustion markers",
    )
    minimal_residual_disease: Optional[float] = Field(
        default=None, ge=0.0, le=1.0,
        description="Estimated MRD fraction from single-cell data",
    )
    actionable_findings: List[str] = Field(
        default_factory=list,
        description="Clinically actionable findings from this timepoint",
    )


# ===================================================================
# PYDANTIC MODELS - WORKFLOW RESULT
# ===================================================================


class WorkflowResult(BaseModel):
    """Container for results from any single-cell analysis workflow."""
    workflow_type: SCWorkflowType = Field(
        ...,
        description="The workflow type that produced these results",
    )
    cell_annotations: List[CellTypeAnnotation] = Field(
        default_factory=list,
        description="Cell type annotation results",
    )
    tme_profile: Optional[TMEProfile] = Field(
        default=None,
        description="Tumor microenvironment profile",
    )
    drug_predictions: List[DrugResponsePrediction] = Field(
        default_factory=list,
        description="Drug response predictions",
    )
    subclones: List[SubclonalResult] = Field(
        default_factory=list,
        description="Subclonal architecture results",
    )
    spatial_niches: List[SpatialNiche] = Field(
        default_factory=list,
        description="Spatial niche results",
    )
    trajectories: List[TrajectoryResult] = Field(
        default_factory=list,
        description="Trajectory analysis results",
    )
    interactions: List[LigandReceptorInteraction] = Field(
        default_factory=list,
        description="Ligand-receptor interaction results",
    )
    biomarkers: List[BiomarkerCandidate] = Field(
        default_factory=list,
        description="Biomarker candidates",
    )
    cart_targets: List[CARTTargetValidation] = Field(
        default_factory=list,
        description="CAR-T target validation results",
    )
    treatment_responses: List[TreatmentResponse] = Field(
        default_factory=list,
        description="Treatment monitoring results",
    )
    severity: SeverityLevel = Field(
        default=SeverityLevel.INFORMATIONAL,
        description="Overall severity of findings",
    )


# ===================================================================
# PYDANTIC MODELS - RESPONSE
# ===================================================================


class SCResponse(BaseModel):
    """Top-level response from the Single-Cell Intelligence Agent."""
    query: str = Field(
        ...,
        description="Original query text",
    )
    workflow_type: SCWorkflowType = Field(
        ...,
        description="Detected or specified workflow type",
    )
    answer: str = Field(
        ...,
        description="LLM-generated answer synthesizing search results",
    )
    search_results: List[SCSearchResult] = Field(
        default_factory=list,
        description="Raw search results from Milvus collections",
    )
    workflow_result: Optional[WorkflowResult] = Field(
        default=None,
        description="Structured workflow-specific results",
    )
    citations: List[str] = Field(
        default_factory=list,
        description="Formatted citation strings",
    )
    confidence: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Overall confidence in the response",
    )
    collections_searched: List[str] = Field(
        default_factory=list,
        description="List of Milvus collections that were queried",
    )
    processing_time_seconds: float = Field(
        default=0.0, ge=0.0,
        description="Total processing time in seconds",
    )


# ===================================================================
# SEARCH PLAN DATACLASS
# ===================================================================


@dataclass
class SearchPlan:
    """Plan for which collections to search and how to weight results.

    Built by the query classifier and consumed by the search orchestrator.
    """
    workflow_type: SCWorkflowType = SCWorkflowType.GENERAL
    collections: List[str] = field(default_factory=list)
    weights: Dict[str, float] = field(default_factory=dict)
    top_k_overrides: Dict[str, int] = field(default_factory=dict)
    filter_expressions: Dict[str, str] = field(default_factory=dict)
    reranking_enabled: bool = True
    explanation: str = ""
