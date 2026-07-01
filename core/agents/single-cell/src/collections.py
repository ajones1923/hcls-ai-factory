"""Milvus collection schemas for the Single-Cell Intelligence Agent.

Defines 12 domain-specific vector collections for single-cell analysis:
  - sc_cell_types       -- Cell type annotations with marker genes
  - sc_markers          -- Gene markers for cell type identification
  - sc_spatial          -- Spatial transcriptomics data and niches
  - sc_tme              -- Tumor microenvironment profiles
  - sc_drug_response    -- Drug sensitivity/resistance predictions
  - sc_literature       -- Published single-cell literature
  - sc_methods          -- Analytical methods and pipelines
  - sc_datasets         -- Reference datasets and atlases
  - sc_trajectories     -- Cellular trajectory and pseudotime data
  - sc_pathways         -- Signaling and metabolic pathways
  - sc_clinical         -- Clinical correlation data
  - genomic_evidence    -- Shared genomic evidence (read-only)

Follows the same pymilvus pattern as:
  rare_disease_diagnostic_agent/src/collections.py

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from pymilvus import (
    CollectionSchema,
    DataType,
    FieldSchema,
)

from src.models import SCWorkflowType


# ===================================================================
# CONSTANTS
# ===================================================================

EMBEDDING_DIM = 384       # BGE-small-en-v1.5
INDEX_TYPE = "IVF_FLAT"
METRIC_TYPE = "COSINE"
NLIST = 128


# ===================================================================
# COLLECTION CONFIG DATACLASS
# ===================================================================


@dataclass
class CollectionConfig:
    """Configuration for a single Milvus vector collection.

    Attributes:
        name: Milvus collection name (e.g. ``sc_cell_types``).
        description: Human-readable description of the collection purpose.
        schema_fields: Ordered list of :class:`pymilvus.FieldSchema` objects
            defining every field in the collection (including id and embedding).
        index_params: Dict of IVF_FLAT / COSINE index parameters.
        estimated_records: Approximate number of records expected after full ingest.
        search_weight: Default relevance weight used during multi-collection search
            (0.0 - 1.0).
    """

    name: str
    description: str
    schema_fields: List[FieldSchema]
    index_params: Dict = field(default_factory=lambda: {
        "metric_type": METRIC_TYPE,
        "index_type": INDEX_TYPE,
        "params": {"nlist": NLIST},
    })
    estimated_records: int = 0
    search_weight: float = 0.05


# ===================================================================
# HELPER -- EMBEDDING FIELD
# ===================================================================


def _make_embedding_field() -> FieldSchema:
    """Create the standard 384-dim FLOAT_VECTOR embedding field.

    All 12 single-cell collections share the same embedding specification
    (BGE-small-en-v1.5, 384 dimensions).

    Returns:
        A :class:`pymilvus.FieldSchema` for the ``embedding`` column.
    """
    return FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding (384-dim)",
    )


# ===================================================================
# COLLECTION SCHEMA DEFINITIONS
# ===================================================================

# -- sc_cell_types -------------------------------------------------

CELL_TYPES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="cell_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Cell type label (e.g., 'CD8+ effector T cell')",
    ),
    FieldSchema(
        name="cell_ontology_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Cell Ontology (CL) identifier",
    ),
    FieldSchema(
        name="tissue",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Tissue of origin",
    ),
    FieldSchema(
        name="compartment",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Cell compartment (immune, stromal, epithelial, etc.)",
    ),
    FieldSchema(
        name="marker_genes",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited canonical marker genes",
    ),
    FieldSchema(
        name="description",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Detailed cell type description and functional role",
    ),
    FieldSchema(
        name="species",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Species (human, mouse)",
    ),
    FieldSchema(
        name="reference_dataset",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Source reference atlas or dataset",
    ),
]

CELL_TYPES_CONFIG = CollectionConfig(
    name="sc_cell_types",
    description="Cell type annotations with ontology mappings, marker genes, and tissue context",
    schema_fields=CELL_TYPES_FIELDS,
    estimated_records=5000,
    search_weight=0.14,
)

# -- sc_markers ----------------------------------------------------

MARKERS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="gene_symbol",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Gene symbol (e.g., CD3E, FOXP3)",
    ),
    FieldSchema(
        name="ensembl_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Ensembl gene identifier",
    ),
    FieldSchema(
        name="cell_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Cell type this marker is associated with",
    ),
    FieldSchema(
        name="tissue",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Tissue context for this marker",
    ),
    FieldSchema(
        name="specificity_score",
        dtype=DataType.FLOAT,
        description="Marker specificity score (0-1, higher = more specific)",
    ),
    FieldSchema(
        name="log2_fold_change",
        dtype=DataType.FLOAT,
        description="Log2 fold change vs other cell types",
    ),
    FieldSchema(
        name="is_surface",
        dtype=DataType.BOOL,
        description="Whether this marker is a cell surface protein",
    ),
    FieldSchema(
        name="evidence_text",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Evidence text supporting this marker association",
    ),
]

MARKERS_CONFIG = CollectionConfig(
    name="sc_markers",
    description="Gene markers for cell type identification with specificity scores and tissue context",
    schema_fields=MARKERS_FIELDS,
    estimated_records=50000,
    search_weight=0.12,
)

# -- sc_spatial ----------------------------------------------------

SPATIAL_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="niche_label",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Spatial niche label (e.g., 'Tumor-immune interface')",
    ),
    FieldSchema(
        name="platform",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Spatial platform (Visium, MERFISH, Xenium, CODEX)",
    ),
    FieldSchema(
        name="tissue",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Tissue type",
    ),
    FieldSchema(
        name="cell_types",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited cell types in this niche",
    ),
    FieldSchema(
        name="signature_genes",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited spatially variable genes",
    ),
    FieldSchema(
        name="morans_i",
        dtype=DataType.FLOAT,
        description="Moran's I spatial autocorrelation statistic",
    ),
    FieldSchema(
        name="clinical_relevance",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Clinical significance of this spatial niche",
    ),
    FieldSchema(
        name="description",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Detailed description of the spatial niche",
    ),
]

SPATIAL_CONFIG = CollectionConfig(
    name="sc_spatial",
    description="Spatial transcriptomics niches with platform-specific data and clinical relevance",
    schema_fields=SPATIAL_FIELDS,
    estimated_records=10000,
    search_weight=0.10,
)

# -- sc_tme --------------------------------------------------------

TME_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="tme_class",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="TME classification (hot_inflamed, cold_desert, excluded, immunosuppressive)",
    ),
    FieldSchema(
        name="cancer_type",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Cancer type or subtype",
    ),
    FieldSchema(
        name="immune_score",
        dtype=DataType.FLOAT,
        description="Aggregate immune infiltration score (0-1)",
    ),
    FieldSchema(
        name="stromal_score",
        dtype=DataType.FLOAT,
        description="Aggregate stromal component score (0-1)",
    ),
    FieldSchema(
        name="cell_composition",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="JSON-encoded cell type composition fractions",
    ),
    FieldSchema(
        name="checkpoint_profile",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited checkpoint gene expression levels",
    ),
    FieldSchema(
        name="immunotherapy_response",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Predicted immunotherapy response category",
    ),
    FieldSchema(
        name="description",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Detailed TME profile description",
    ),
]

TME_CONFIG = CollectionConfig(
    name="sc_tme",
    description="Tumor microenvironment profiles with immune phenotyping and therapy prediction",
    schema_fields=TME_FIELDS,
    estimated_records=8000,
    search_weight=0.10,
)

# -- sc_drug_response ----------------------------------------------

DRUG_RESPONSE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="drug_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Drug or compound name",
    ),
    FieldSchema(
        name="drug_class",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Drug class (e.g., checkpoint inhibitor, TKI)",
    ),
    FieldSchema(
        name="cancer_type",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Cancer type evaluated",
    ),
    FieldSchema(
        name="cell_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Target cell type or subpopulation",
    ),
    FieldSchema(
        name="sensitivity_score",
        dtype=DataType.FLOAT,
        description="Predicted sensitivity score (0-1)",
    ),
    FieldSchema(
        name="resistance_mechanisms",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited resistance mechanisms",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Data source (DepMap, CCLE, GDSC, clinical)",
    ),
    FieldSchema(
        name="evidence_text",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Supporting evidence text",
    ),
]

DRUG_RESPONSE_CONFIG = CollectionConfig(
    name="sc_drug_response",
    description="Drug sensitivity/resistance predictions from single-cell transcriptomic signatures",
    schema_fields=DRUG_RESPONSE_FIELDS,
    estimated_records=25000,
    search_weight=0.09,
)

# -- sc_literature -------------------------------------------------

LITERATURE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="pmid",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="PubMed identifier",
    ),
    FieldSchema(
        name="doi",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Digital Object Identifier",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Publication title",
    ),
    FieldSchema(
        name="authors",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Author list (pipe-delimited)",
    ),
    FieldSchema(
        name="journal",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Journal name",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Publication year",
    ),
    FieldSchema(
        name="abstract_text",
        dtype=DataType.VARCHAR,
        max_length=8192,
        description="Publication abstract text",
    ),
    FieldSchema(
        name="topics",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Pipe-delimited topic tags",
    ),
]

LITERATURE_CONFIG = CollectionConfig(
    name="sc_literature",
    description="Published single-cell transcriptomics literature with abstracts and topic tags",
    schema_fields=LITERATURE_FIELDS,
    estimated_records=50000,
    search_weight=0.08,
)

# -- sc_methods ----------------------------------------------------

METHODS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="method_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Method or tool name (e.g., Scanpy, Seurat, CellChat)",
    ),
    FieldSchema(
        name="category",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Method category (preprocessing, clustering, trajectory, etc.)",
    ),
    FieldSchema(
        name="assay_types",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Pipe-delimited compatible assay types",
    ),
    FieldSchema(
        name="description",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Method description and use cases",
    ),
    FieldSchema(
        name="strengths",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Method strengths and best-suited scenarios",
    ),
    FieldSchema(
        name="limitations",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Known limitations and caveats",
    ),
    FieldSchema(
        name="reference_pmid",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="PubMed ID of the method publication",
    ),
    FieldSchema(
        name="gpu_accelerated",
        dtype=DataType.BOOL,
        description="Whether the method supports GPU acceleration (e.g., RAPIDS)",
    ),
]

METHODS_CONFIG = CollectionConfig(
    name="sc_methods",
    description="Analytical methods and computational tools for single-cell analysis",
    schema_fields=METHODS_FIELDS,
    estimated_records=2000,
    search_weight=0.07,
)

# -- sc_datasets ---------------------------------------------------

DATASETS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="dataset_id",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Dataset identifier (CellxGene, GEO accession, etc.)",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Dataset name or title",
    ),
    FieldSchema(
        name="tissue",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Tissue type",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Disease context (or 'normal')",
    ),
    FieldSchema(
        name="assay_type",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Assay modality (scRNA-seq, CITE-seq, etc.)",
    ),
    FieldSchema(
        name="cell_count",
        dtype=DataType.INT64,
        description="Number of cells in the dataset",
    ),
    FieldSchema(
        name="cell_types",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Pipe-delimited cell types present in the dataset",
    ),
    FieldSchema(
        name="species",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Species (human, mouse)",
    ),
    FieldSchema(
        name="description",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Dataset description and experimental design",
    ),
]

DATASETS_CONFIG = CollectionConfig(
    name="sc_datasets",
    description="Reference single-cell datasets and atlases with metadata and cell type inventories",
    schema_fields=DATASETS_FIELDS,
    estimated_records=15000,
    search_weight=0.06,
)

# -- sc_trajectories -----------------------------------------------

TRAJECTORIES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="trajectory_type",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Trajectory type (differentiation, activation, exhaustion, EMT, etc.)",
    ),
    FieldSchema(
        name="start_cell_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Starting cell state",
    ),
    FieldSchema(
        name="end_cell_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Terminal cell state",
    ),
    FieldSchema(
        name="tissue",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Tissue context",
    ),
    FieldSchema(
        name="driver_genes",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited genes driving the trajectory",
    ),
    FieldSchema(
        name="method",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Trajectory inference method (Monocle, PAGA, RNA velocity, etc.)",
    ),
    FieldSchema(
        name="clinical_relevance",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Clinical significance of the trajectory",
    ),
    FieldSchema(
        name="description",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Detailed trajectory description",
    ),
]

TRAJECTORIES_CONFIG = CollectionConfig(
    name="sc_trajectories",
    description="Cellular trajectory and pseudotime analysis results with driver genes",
    schema_fields=TRAJECTORIES_FIELDS,
    estimated_records=8000,
    search_weight=0.07,
)

# -- sc_pathways ---------------------------------------------------

PATHWAYS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="pathway_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Pathway identifier (KEGG, Reactome, GO)",
    ),
    FieldSchema(
        name="pathway_name",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Pathway name",
    ),
    FieldSchema(
        name="source_db",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Source database (KEGG, Reactome, MSigDB, GO)",
    ),
    FieldSchema(
        name="genes",
        dtype=DataType.VARCHAR,
        max_length=8192,
        description="Pipe-delimited pathway member genes",
    ),
    FieldSchema(
        name="cell_type_activity",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="JSON-encoded cell-type-specific pathway activity scores",
    ),
    FieldSchema(
        name="description",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Pathway description and biological function",
    ),
]

PATHWAYS_CONFIG = CollectionConfig(
    name="sc_pathways",
    description="Signaling and metabolic pathways with cell-type-specific activity profiles",
    schema_fields=PATHWAYS_FIELDS,
    estimated_records=20000,
    search_weight=0.07,
)

# -- sc_clinical ---------------------------------------------------

CLINICAL_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="indication",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Clinical indication or disease",
    ),
    FieldSchema(
        name="biomarker_gene",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Biomarker gene symbol",
    ),
    FieldSchema(
        name="biomarker_type",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Biomarker type (diagnostic, prognostic, predictive, pharmacodynamic)",
    ),
    FieldSchema(
        name="cell_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Relevant cell type",
    ),
    FieldSchema(
        name="clinical_outcome",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Associated clinical outcome",
    ),
    FieldSchema(
        name="hazard_ratio",
        dtype=DataType.FLOAT,
        description="Hazard ratio for prognostic markers (0 if not applicable)",
    ),
    FieldSchema(
        name="trial_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Associated clinical trial NCT ID",
    ),
    FieldSchema(
        name="evidence_text",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Clinical evidence summary",
    ),
]

CLINICAL_CONFIG = CollectionConfig(
    name="sc_clinical",
    description="Clinical correlation data linking single-cell biomarkers to patient outcomes",
    schema_fields=CLINICAL_FIELDS,
    estimated_records=12000,
    search_weight=0.07,
)

# -- genomic_evidence (shared, read-only) --------------------------

GENOMIC_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Gene symbol",
    ),
    FieldSchema(
        name="variant",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Variant notation (HGVS)",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Evidence source (ClinVar, AlphaMissense, etc.)",
    ),
    FieldSchema(
        name="classification",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Variant classification",
    ),
    FieldSchema(
        name="evidence_text",
        dtype=DataType.VARCHAR,
        max_length=8192,
        description="Evidence text from ClinVar, AlphaMissense, or literature",
    ),
]

GENOMIC_CONFIG = CollectionConfig(
    name="genomic_evidence",
    description="Shared genomic evidence collection (read-only, managed by genomics pipeline)",
    schema_fields=GENOMIC_FIELDS,
    estimated_records=3560000,
    search_weight=0.03,
)


# ===================================================================
# ALL COLLECTIONS LIST
# ===================================================================

ALL_COLLECTIONS: List[CollectionConfig] = [
    CELL_TYPES_CONFIG,
    MARKERS_CONFIG,
    SPATIAL_CONFIG,
    TME_CONFIG,
    DRUG_RESPONSE_CONFIG,
    LITERATURE_CONFIG,
    METHODS_CONFIG,
    DATASETS_CONFIG,
    TRAJECTORIES_CONFIG,
    PATHWAYS_CONFIG,
    CLINICAL_CONFIG,
    GENOMIC_CONFIG,
]
"""Ordered list of all 12 single-cell collection configurations."""


COLLECTION_NAMES: Dict[str, str] = {
    "cell_types": "sc_cell_types",
    "markers": "sc_markers",
    "spatial": "sc_spatial",
    "tme": "sc_tme",
    "drug_response": "sc_drug_response",
    "literature": "sc_literature",
    "methods": "sc_methods",
    "datasets": "sc_datasets",
    "trajectories": "sc_trajectories",
    "pathways": "sc_pathways",
    "clinical": "sc_clinical",
    "genomic": "genomic_evidence",
}
"""Short-name to full-name mapping for convenience lookups."""

_FULL_NAME_MAP: Dict[str, CollectionConfig] = {
    cfg.name: cfg for cfg in ALL_COLLECTIONS
}


# ===================================================================
# DEFAULT SEARCH WEIGHTS
# ===================================================================

_DEFAULT_SEARCH_WEIGHTS: Dict[str, float] = {
    cfg.name: cfg.search_weight for cfg in ALL_COLLECTIONS
}
"""Base search weights used when no workflow-specific boost is applied.
Sum: {sum:.2f}.""".format(sum=sum(cfg.search_weight for cfg in ALL_COLLECTIONS))


# ===================================================================
# WORKFLOW COLLECTION WEIGHTS
# ===================================================================

WORKFLOW_COLLECTION_WEIGHTS: Dict[SCWorkflowType, Dict[str, float]] = {

    # -- Cell Type Annotation --------------------------------------
    SCWorkflowType.CELL_TYPE_ANNOTATION: {
        "sc_cell_types": 0.25,
        "sc_markers": 0.22,
        "sc_datasets": 0.10,
        "sc_methods": 0.10,
        "sc_literature": 0.08,
        "sc_pathways": 0.06,
        "sc_trajectories": 0.05,
        "sc_tme": 0.04,
        "sc_spatial": 0.03,
        "sc_clinical": 0.03,
        "sc_drug_response": 0.02,
        "genomic_evidence": 0.02,
    },

    # -- TME Profiling ---------------------------------------------
    SCWorkflowType.TME_PROFILING: {
        "sc_tme": 0.25,
        "sc_cell_types": 0.15,
        "sc_markers": 0.10,
        "sc_spatial": 0.10,
        "sc_drug_response": 0.08,
        "sc_literature": 0.08,
        "sc_clinical": 0.06,
        "sc_pathways": 0.05,
        "sc_trajectories": 0.04,
        "sc_datasets": 0.04,
        "sc_methods": 0.03,
        "genomic_evidence": 0.02,
    },

    # -- Drug Response ---------------------------------------------
    SCWorkflowType.DRUG_RESPONSE: {
        "sc_drug_response": 0.25,
        "sc_tme": 0.12,
        "sc_cell_types": 0.10,
        "sc_clinical": 0.10,
        "sc_markers": 0.08,
        "sc_literature": 0.08,
        "sc_pathways": 0.07,
        "sc_trajectories": 0.05,
        "sc_spatial": 0.04,
        "sc_datasets": 0.04,
        "sc_methods": 0.04,
        "genomic_evidence": 0.03,
    },

    # -- Subclonal Architecture ------------------------------------
    SCWorkflowType.SUBCLONAL_ARCHITECTURE: {
        "sc_cell_types": 0.15,
        "sc_markers": 0.10,
        "genomic_evidence": 0.15,
        "sc_drug_response": 0.12,
        "sc_trajectories": 0.10,
        "sc_tme": 0.08,
        "sc_literature": 0.08,
        "sc_clinical": 0.06,
        "sc_pathways": 0.05,
        "sc_methods": 0.04,
        "sc_spatial": 0.04,
        "sc_datasets": 0.03,
    },

    # -- Spatial Niche ---------------------------------------------
    SCWorkflowType.SPATIAL_NICHE: {
        "sc_spatial": 0.28,
        "sc_cell_types": 0.12,
        "sc_tme": 0.12,
        "sc_markers": 0.08,
        "sc_methods": 0.08,
        "sc_literature": 0.08,
        "sc_pathways": 0.06,
        "sc_clinical": 0.05,
        "sc_drug_response": 0.04,
        "sc_datasets": 0.04,
        "sc_trajectories": 0.03,
        "genomic_evidence": 0.02,
    },

    # -- Trajectory Analysis ---------------------------------------
    SCWorkflowType.TRAJECTORY_ANALYSIS: {
        "sc_trajectories": 0.25,
        "sc_cell_types": 0.15,
        "sc_markers": 0.10,
        "sc_methods": 0.10,
        "sc_pathways": 0.08,
        "sc_literature": 0.08,
        "sc_datasets": 0.06,
        "sc_tme": 0.05,
        "sc_spatial": 0.04,
        "sc_clinical": 0.04,
        "sc_drug_response": 0.03,
        "genomic_evidence": 0.02,
    },

    # -- Ligand-Receptor Interactions ------------------------------
    SCWorkflowType.LIGAND_RECEPTOR: {
        "sc_cell_types": 0.15,
        "sc_markers": 0.15,
        "sc_pathways": 0.15,
        "sc_tme": 0.10,
        "sc_spatial": 0.10,
        "sc_methods": 0.08,
        "sc_literature": 0.08,
        "sc_drug_response": 0.05,
        "sc_clinical": 0.04,
        "sc_trajectories": 0.04,
        "sc_datasets": 0.04,
        "genomic_evidence": 0.02,
    },

    # -- Biomarker Discovery ---------------------------------------
    SCWorkflowType.BIOMARKER_DISCOVERY: {
        "sc_markers": 0.20,
        "sc_clinical": 0.18,
        "sc_cell_types": 0.12,
        "sc_literature": 0.10,
        "sc_drug_response": 0.08,
        "sc_tme": 0.06,
        "sc_pathways": 0.06,
        "sc_datasets": 0.06,
        "sc_trajectories": 0.04,
        "sc_spatial": 0.04,
        "sc_methods": 0.03,
        "genomic_evidence": 0.03,
    },

    # -- CAR-T Target Validation -----------------------------------
    SCWorkflowType.CART_TARGET_VALIDATION: {
        "sc_markers": 0.18,
        "sc_cell_types": 0.15,
        "sc_tme": 0.12,
        "sc_clinical": 0.12,
        "sc_drug_response": 0.10,
        "sc_literature": 0.08,
        "sc_datasets": 0.06,
        "sc_spatial": 0.05,
        "sc_pathways": 0.05,
        "genomic_evidence": 0.04,
        "sc_trajectories": 0.03,
        "sc_methods": 0.02,
    },

    # -- Treatment Monitoring --------------------------------------
    SCWorkflowType.TREATMENT_MONITORING: {
        "sc_clinical": 0.20,
        "sc_drug_response": 0.18,
        "sc_cell_types": 0.12,
        "sc_tme": 0.10,
        "sc_markers": 0.08,
        "sc_trajectories": 0.08,
        "sc_literature": 0.06,
        "sc_pathways": 0.05,
        "sc_spatial": 0.04,
        "sc_datasets": 0.04,
        "sc_methods": 0.03,
        "genomic_evidence": 0.02,
    },

    # -- General ---------------------------------------------------
    SCWorkflowType.GENERAL: {
        "sc_cell_types": 0.14,
        "sc_markers": 0.12,
        "sc_spatial": 0.10,
        "sc_tme": 0.10,
        "sc_drug_response": 0.09,
        "sc_literature": 0.08,
        "sc_methods": 0.07,
        "sc_datasets": 0.06,
        "sc_trajectories": 0.07,
        "sc_pathways": 0.07,
        "sc_clinical": 0.07,
        "genomic_evidence": 0.03,
    },
}


# ===================================================================
# HELPER FUNCTIONS
# ===================================================================


def get_collection_config(name: str) -> CollectionConfig:
    """Look up a :class:`CollectionConfig` by full collection name.

    Args:
        name: Full Milvus collection name (e.g. ``sc_cell_types``)
            **or** a short alias (e.g. ``cell_types``).

    Returns:
        The matching :class:`CollectionConfig`.

    Raises:
        KeyError: If no collection matches *name*.
    """
    # Try full name first
    if name in _FULL_NAME_MAP:
        return _FULL_NAME_MAP[name]

    # Try short alias
    full_name = COLLECTION_NAMES.get(name)
    if full_name and full_name in _FULL_NAME_MAP:
        return _FULL_NAME_MAP[full_name]

    valid = ", ".join(sorted(_FULL_NAME_MAP.keys()))
    raise KeyError(
        f"Unknown collection '{name}'. "
        f"Valid collections: {valid}"
    )


def get_all_collection_names() -> List[str]:
    """Return a list of all 12 full Milvus collection names.

    Returns:
        Ordered list of collection name strings.
    """
    return [cfg.name for cfg in ALL_COLLECTIONS]


def get_search_weights(
    workflow: Optional[SCWorkflowType] = None,
) -> Dict[str, float]:
    """Return collection search weights, optionally boosted for a workflow.

    When *workflow* is ``None`` (or not found in the boost table), the
    default weights from each :class:`CollectionConfig` are returned.

    Args:
        workflow: Optional workflow type to apply workflow-specific weighting.

    Returns:
        Dict mapping collection name to its search weight (0.0 - 1.0).
    """
    if workflow is not None and workflow in WORKFLOW_COLLECTION_WEIGHTS:
        return dict(WORKFLOW_COLLECTION_WEIGHTS[workflow])
    return dict(_DEFAULT_SEARCH_WEIGHTS)


def get_collection_schema(name: str) -> CollectionSchema:
    """Build a :class:`pymilvus.CollectionSchema` for the named collection.

    Args:
        name: Full or short collection name.

    Returns:
        A ready-to-use :class:`CollectionSchema`.
    """
    cfg = get_collection_config(name)
    return CollectionSchema(
        fields=cfg.schema_fields,
        description=cfg.description,
    )
