"""Milvus collection schemas for the Rare Disease Diagnostic Agent.

Defines 14 domain-specific vector collections for rare disease diagnostics:
  - rd_phenotypes          -- HPO phenotype terms with information content
  - rd_diseases            -- Rare disease entries (OMIM/Orphanet)
  - rd_genes               -- Disease-associated genes with constraint scores
  - rd_variants            -- Classified genetic variants (ACMG)
  - rd_literature          -- Published rare disease literature
  - rd_trials              -- Clinical trials for rare disease therapies
  - rd_therapies           -- Approved and investigational therapies
  - rd_case_reports        -- Published case reports with phenotype-genotype data
  - rd_guidelines          -- Clinical practice guidelines
  - rd_pathways            -- Metabolic and signaling pathways
  - rd_registries          -- Patient registries and natural history studies
  - rd_natural_history     -- Disease natural history milestones
  - rd_newborn_screening   -- Newborn screening conditions and protocols
  - genomic_evidence       -- Shared genomic evidence (read-only)

Follows the same pymilvus pattern as:
  clinical_trial_intelligence_agent/src/collections.py

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

from src.models import DiagnosticWorkflowType


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
        name: Milvus collection name (e.g. ``rd_phenotypes``).
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

    All 14 rare disease collections share the same embedding specification
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

# -- rd_phenotypes ------------------------------------------------

PHENOTYPES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="hpo_id",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="HPO term identifier (e.g., HP:0001250)",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="HPO term name (e.g., Seizure)",
    ),
    FieldSchema(
        name="definition",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="HPO term definition text",
    ),
    FieldSchema(
        name="synonyms",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited synonyms for the phenotype",
    ),
    FieldSchema(
        name="ic_score",
        dtype=DataType.FLOAT,
        description="Information content score (higher = more specific)",
    ),
    FieldSchema(
        name="frequency",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Phenotype frequency class (e.g., HP:0040281 Very frequent)",
    ),
    FieldSchema(
        name="is_negated",
        dtype=DataType.BOOL,
        description="Whether this is a negated (absent) phenotype observation",
    ),
]

PHENOTYPES_CONFIG = CollectionConfig(
    name="rd_phenotypes",
    description="HPO phenotype terms with definitions, synonyms, and information content scores",
    schema_fields=PHENOTYPES_FIELDS,
    estimated_records=18000,
    search_weight=0.12,
)

# -- rd_diseases --------------------------------------------------

DISEASES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="disease_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Disease identifier (OMIM or ORPHA number)",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Disease name",
    ),
    FieldSchema(
        name="omim_id",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="OMIM disease number",
    ),
    FieldSchema(
        name="orpha_code",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Orphanet disease code",
    ),
    FieldSchema(
        name="inheritance",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Inheritance pattern (AD, AR, XL, MT, etc.)",
    ),
    FieldSchema(
        name="prevalence",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Prevalence estimate (e.g., 1/100000)",
    ),
    FieldSchema(
        name="category",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Disease category (metabolic, neurological, etc.)",
    ),
    FieldSchema(
        name="clinical_features",
        dtype=DataType.VARCHAR,
        max_length=8192,
        description="Key clinical features and HPO associations as text",
    ),
]

DISEASES_CONFIG = CollectionConfig(
    name="rd_diseases",
    description="Rare disease entries from OMIM and Orphanet with inheritance and clinical features",
    schema_fields=DISEASES_FIELDS,
    estimated_records=10000,
    search_weight=0.11,
)

# -- rd_genes -----------------------------------------------------

GENES_FIELDS = [
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
        description="HGNC gene symbol (e.g., CFTR)",
    ),
    FieldSchema(
        name="gene_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Full gene name",
    ),
    FieldSchema(
        name="chromosome",
        dtype=DataType.VARCHAR,
        max_length=8,
        description="Chromosome location (e.g., 7q31.2)",
    ),
    FieldSchema(
        name="constraint_score",
        dtype=DataType.FLOAT,
        description="gnomAD constraint score (pLI or LOEUF)",
    ),
    FieldSchema(
        name="disease_associations",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Pipe-delimited disease associations",
    ),
    FieldSchema(
        name="function_summary",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Gene function summary from UniProt/NCBI",
    ),
]

GENES_CONFIG = CollectionConfig(
    name="rd_genes",
    description="Disease-associated genes with constraint scores and functional summaries",
    schema_fields=GENES_FIELDS,
    estimated_records=5000,
    search_weight=0.10,
)

# -- rd_variants --------------------------------------------------

VARIANTS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="variant_id",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Variant identifier (chr:pos:ref:alt or dbSNP rsID)",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Gene symbol (HGNC)",
    ),
    FieldSchema(
        name="hgvs",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="HGVS notation (coding or protein)",
    ),
    FieldSchema(
        name="classification",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ACMG classification (pathogenic, likely_pathogenic, vus, etc.)",
    ),
    FieldSchema(
        name="population_freq",
        dtype=DataType.FLOAT,
        description="Population allele frequency from gnomAD",
    ),
    FieldSchema(
        name="clinvar_stars",
        dtype=DataType.INT32,
        description="ClinVar review star rating (0-4)",
    ),
    FieldSchema(
        name="review_status",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="ClinVar review status text",
    ),
]

VARIANTS_CONFIG = CollectionConfig(
    name="rd_variants",
    description="ACMG-classified genetic variants with ClinVar review status and population frequencies",
    schema_fields=VARIANTS_FIELDS,
    estimated_records=500000,
    search_weight=0.10,
)

# -- rd_literature ------------------------------------------------

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
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Article title",
    ),
    FieldSchema(
        name="abstract",
        dtype=DataType.VARCHAR,
        max_length=8192,
        description="Article abstract text",
    ),
    FieldSchema(
        name="journal",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Journal name",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT32,
        description="Publication year",
    ),
    FieldSchema(
        name="disease_context",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Primary disease context of the article",
    ),
]

LITERATURE_CONFIG = CollectionConfig(
    name="rd_literature",
    description="Published rare disease literature with abstracts and disease context",
    schema_fields=LITERATURE_FIELDS,
    estimated_records=50000,
    search_weight=0.08,
)

# -- rd_trials ----------------------------------------------------

TRIALS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="nct_id",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Trial title",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Target condition(s)",
    ),
    FieldSchema(
        name="intervention",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Intervention description",
    ),
    FieldSchema(
        name="phase",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Trial phase (Phase 1, Phase 2, Phase 3, etc.)",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Trial recruitment status",
    ),
    FieldSchema(
        name="eligibility",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Key eligibility criteria text",
    ),
]

TRIALS_CONFIG = CollectionConfig(
    name="rd_trials",
    description="Clinical trials targeting rare diseases with eligibility criteria",
    schema_fields=TRIALS_FIELDS,
    estimated_records=8000,
    search_weight=0.06,
)

# -- rd_therapies -------------------------------------------------

THERAPIES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="therapy_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Drug or therapy name",
    ),
    FieldSchema(
        name="indication",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Approved indication",
    ),
    FieldSchema(
        name="mechanism",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Mechanism of action",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Approval status (approved_fda, approved_ema, investigational, etc.)",
    ),
    FieldSchema(
        name="approval_year",
        dtype=DataType.INT32,
        description="Year of regulatory approval (0 if not yet approved)",
    ),
    FieldSchema(
        name="gene_target",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Target gene for gene therapies",
    ),
]

THERAPIES_CONFIG = CollectionConfig(
    name="rd_therapies",
    description="Approved and investigational rare disease therapies including gene therapies",
    schema_fields=THERAPIES_FIELDS,
    estimated_records=2000,
    search_weight=0.07,
)

# -- rd_case_reports ----------------------------------------------

CASE_REPORTS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="case_id",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Case report identifier",
    ),
    FieldSchema(
        name="phenotypes_hpo",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Pipe-delimited HPO term IDs observed in the case",
    ),
    FieldSchema(
        name="diagnosis",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Final diagnosis",
    ),
    FieldSchema(
        name="genotype",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Genotype findings (gene:variant notation)",
    ),
    FieldSchema(
        name="age_onset",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Age at symptom onset",
    ),
    FieldSchema(
        name="outcome",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Clinical outcome and follow-up",
    ),
]

CASE_REPORTS_CONFIG = CollectionConfig(
    name="rd_case_reports",
    description="Published rare disease case reports with phenotype, genotype, and outcome data",
    schema_fields=CASE_REPORTS_FIELDS,
    estimated_records=20000,
    search_weight=0.07,
)

# -- rd_guidelines ------------------------------------------------

GUIDELINES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="guideline_id",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Guideline identifier",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Guideline title",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Publishing organization (ACMG, ESHG, GeneReviews, etc.)",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Target disease or condition",
    ),
    FieldSchema(
        name="recommendation",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Recommendation text chunk",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Evidence level (strong, moderate, limited, etc.)",
    ),
]

GUIDELINES_CONFIG = CollectionConfig(
    name="rd_guidelines",
    description="Clinical practice guidelines for rare disease diagnosis and management",
    schema_fields=GUIDELINES_FIELDS,
    estimated_records=3000,
    search_weight=0.06,
)

# -- rd_pathways --------------------------------------------------

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
        description="Pathway identifier (KEGG, Reactome, etc.)",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Pathway name",
    ),
    FieldSchema(
        name="genes",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Pipe-delimited gene symbols in the pathway",
    ),
    FieldSchema(
        name="enzymes",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited enzyme names/EC numbers",
    ),
    FieldSchema(
        name="metabolites",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Pipe-delimited metabolite names",
    ),
    FieldSchema(
        name="disease_associations",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Pipe-delimited associated diseases",
    ),
]

PATHWAYS_CONFIG = CollectionConfig(
    name="rd_pathways",
    description="Metabolic and signaling pathways with gene, enzyme, and disease associations",
    schema_fields=PATHWAYS_FIELDS,
    estimated_records=2000,
    search_weight=0.06,
)

# -- rd_registries ------------------------------------------------

REGISTRIES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="registry_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Patient registry name",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Target disease or disease group",
    ),
    FieldSchema(
        name="organization",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Sponsoring organization",
    ),
    FieldSchema(
        name="enrollment",
        dtype=DataType.INT32,
        description="Current enrollment count",
    ),
    FieldSchema(
        name="country",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Primary country or 'international'",
    ),
]

REGISTRIES_CONFIG = CollectionConfig(
    name="rd_registries",
    description="Patient registries for rare diseases with enrollment and organization data",
    schema_fields=REGISTRIES_FIELDS,
    estimated_records=1500,
    search_weight=0.04,
)

# -- rd_natural_history -------------------------------------------

NATURAL_HISTORY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Disease name",
    ),
    FieldSchema(
        name="milestone",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Natural history milestone or event",
    ),
    FieldSchema(
        name="age_range",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Typical age range for the milestone",
    ),
    FieldSchema(
        name="frequency",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Frequency of the milestone in affected patients",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Data source (registry, cohort study, literature review)",
    ),
]

NATURAL_HISTORY_CONFIG = CollectionConfig(
    name="rd_natural_history",
    description="Disease natural history milestones with age ranges and frequencies",
    schema_fields=NATURAL_HISTORY_FIELDS,
    estimated_records=5000,
    search_weight=0.05,
)

# -- rd_newborn_screening -----------------------------------------

NEWBORN_SCREENING_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Screened condition name",
    ),
    FieldSchema(
        name="analyte",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Screening analyte or biomarker",
    ),
    FieldSchema(
        name="cutoff",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Screening cutoff value and units",
    ),
    FieldSchema(
        name="confirmatory_test",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Recommended confirmatory testing",
    ),
    FieldSchema(
        name="act_sheet",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="ACT sheet content (ACMG newborn screening ACT sheets)",
    ),
]

NEWBORN_SCREENING_CONFIG = CollectionConfig(
    name="rd_newborn_screening",
    description="Newborn screening conditions with analytes, cutoffs, and ACMG ACT sheets",
    schema_fields=NEWBORN_SCREENING_FIELDS,
    estimated_records=80,
    search_weight=0.05,
)

# -- genomic_evidence (shared, read-only) -------------------------

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
        name="variant_id",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Genomic variant identifier",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Gene symbol",
    ),
    FieldSchema(
        name="consequence",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Variant consequence (missense, nonsense, frameshift, etc.)",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="ClinVar clinical significance",
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
    PHENOTYPES_CONFIG,
    DISEASES_CONFIG,
    GENES_CONFIG,
    VARIANTS_CONFIG,
    LITERATURE_CONFIG,
    TRIALS_CONFIG,
    THERAPIES_CONFIG,
    CASE_REPORTS_CONFIG,
    GUIDELINES_CONFIG,
    PATHWAYS_CONFIG,
    REGISTRIES_CONFIG,
    NATURAL_HISTORY_CONFIG,
    NEWBORN_SCREENING_CONFIG,
    GENOMIC_CONFIG,
]
"""Ordered list of all 14 rare disease collection configurations."""


COLLECTION_NAMES: Dict[str, str] = {
    "phenotypes": "rd_phenotypes",
    "diseases": "rd_diseases",
    "genes": "rd_genes",
    "variants": "rd_variants",
    "literature": "rd_literature",
    "trials": "rd_trials",
    "therapies": "rd_therapies",
    "case_reports": "rd_case_reports",
    "guidelines": "rd_guidelines",
    "pathways": "rd_pathways",
    "registries": "rd_registries",
    "natural_history": "rd_natural_history",
    "newborn_screening": "rd_newborn_screening",
    "genomic": "genomic_evidence",
}
"""Mapping of short alias names to full Milvus collection names."""


# ===================================================================
# COLLECTION SCHEMAS (pymilvus CollectionSchema objects)
# ===================================================================

COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    cfg.name: CollectionSchema(
        fields=cfg.schema_fields,
        description=cfg.description,
    )
    for cfg in ALL_COLLECTIONS
}
"""Mapping of collection name to pymilvus CollectionSchema, ready for
``Collection(name=..., schema=...)`` instantiation."""


# ===================================================================
# DEFAULT SEARCH WEIGHTS
# ===================================================================

_DEFAULT_SEARCH_WEIGHTS: Dict[str, float] = {
    cfg.name: cfg.search_weight for cfg in ALL_COLLECTIONS
}
"""Base search weights used when no workflow-specific boost is applied.
Sum: {sum:.2f}.""".format(sum=sum(cfg.search_weight for cfg in ALL_COLLECTIONS))


# ===================================================================
# WORKFLOW-SPECIFIC COLLECTION WEIGHTS
# ===================================================================

WORKFLOW_COLLECTION_WEIGHTS: Dict[DiagnosticWorkflowType, Dict[str, float]] = {

    # -- Phenotype-Driven Diagnosis --------------------------------
    DiagnosticWorkflowType.PHENOTYPE_DRIVEN: {
        "rd_phenotypes": 0.22,
        "rd_diseases": 0.18,
        "rd_case_reports": 0.12,
        "rd_genes": 0.10,
        "rd_literature": 0.08,
        "rd_natural_history": 0.06,
        "rd_guidelines": 0.05,
        "rd_variants": 0.05,
        "rd_pathways": 0.04,
        "rd_therapies": 0.03,
        "rd_registries": 0.02,
        "rd_trials": 0.02,
        "rd_newborn_screening": 0.02,
        "genomic_evidence": 0.01,
    },

    # -- WES/WGS Interpretation ------------------------------------
    DiagnosticWorkflowType.WES_WGS_INTERPRETATION: {
        "rd_variants": 0.22,
        "rd_genes": 0.18,
        "rd_diseases": 0.12,
        "rd_literature": 0.10,
        "rd_phenotypes": 0.08,
        "rd_guidelines": 0.06,
        "rd_case_reports": 0.05,
        "genomic_evidence": 0.05,
        "rd_pathways": 0.04,
        "rd_therapies": 0.03,
        "rd_natural_history": 0.03,
        "rd_trials": 0.02,
        "rd_registries": 0.01,
        "rd_newborn_screening": 0.01,
    },

    # -- Metabolic Screening ---------------------------------------
    DiagnosticWorkflowType.METABOLIC_SCREENING: {
        "rd_pathways": 0.20,
        "rd_newborn_screening": 0.15,
        "rd_diseases": 0.12,
        "rd_genes": 0.10,
        "rd_phenotypes": 0.08,
        "rd_guidelines": 0.07,
        "rd_literature": 0.06,
        "rd_case_reports": 0.05,
        "rd_variants": 0.05,
        "rd_therapies": 0.04,
        "rd_natural_history": 0.03,
        "rd_registries": 0.02,
        "rd_trials": 0.02,
        "genomic_evidence": 0.01,
    },

    # -- Dysmorphology ---------------------------------------------
    DiagnosticWorkflowType.DYSMORPHOLOGY: {
        "rd_phenotypes": 0.25,
        "rd_diseases": 0.18,
        "rd_case_reports": 0.12,
        "rd_genes": 0.10,
        "rd_literature": 0.08,
        "rd_guidelines": 0.06,
        "rd_natural_history": 0.05,
        "rd_variants": 0.04,
        "rd_pathways": 0.03,
        "rd_registries": 0.03,
        "rd_therapies": 0.02,
        "rd_trials": 0.02,
        "rd_newborn_screening": 0.01,
        "genomic_evidence": 0.01,
    },

    # -- Neurogenetic Evaluation -----------------------------------
    DiagnosticWorkflowType.NEUROGENETIC: {
        "rd_genes": 0.18,
        "rd_diseases": 0.15,
        "rd_phenotypes": 0.12,
        "rd_variants": 0.12,
        "rd_literature": 0.08,
        "rd_case_reports": 0.08,
        "rd_natural_history": 0.06,
        "rd_guidelines": 0.05,
        "rd_pathways": 0.05,
        "rd_therapies": 0.04,
        "rd_trials": 0.03,
        "genomic_evidence": 0.02,
        "rd_registries": 0.01,
        "rd_newborn_screening": 0.01,
    },

    # -- Cardiac Genetics ------------------------------------------
    DiagnosticWorkflowType.CARDIAC_GENETICS: {
        "rd_genes": 0.18,
        "rd_variants": 0.15,
        "rd_diseases": 0.14,
        "rd_phenotypes": 0.10,
        "rd_guidelines": 0.10,
        "rd_literature": 0.08,
        "rd_case_reports": 0.06,
        "rd_natural_history": 0.05,
        "rd_therapies": 0.04,
        "rd_pathways": 0.03,
        "rd_trials": 0.03,
        "genomic_evidence": 0.02,
        "rd_registries": 0.01,
        "rd_newborn_screening": 0.01,
    },

    # -- Connective Tissue Disorders -------------------------------
    DiagnosticWorkflowType.CONNECTIVE_TISSUE: {
        "rd_phenotypes": 0.18,
        "rd_diseases": 0.15,
        "rd_genes": 0.14,
        "rd_guidelines": 0.10,
        "rd_variants": 0.08,
        "rd_case_reports": 0.08,
        "rd_literature": 0.07,
        "rd_natural_history": 0.06,
        "rd_pathways": 0.04,
        "rd_therapies": 0.03,
        "rd_trials": 0.03,
        "rd_registries": 0.02,
        "genomic_evidence": 0.01,
        "rd_newborn_screening": 0.01,
    },

    # -- Inborn Errors of Metabolism -------------------------------
    DiagnosticWorkflowType.INBORN_ERRORS: {
        "rd_pathways": 0.20,
        "rd_diseases": 0.14,
        "rd_genes": 0.12,
        "rd_newborn_screening": 0.10,
        "rd_phenotypes": 0.08,
        "rd_guidelines": 0.08,
        "rd_variants": 0.06,
        "rd_literature": 0.05,
        "rd_therapies": 0.05,
        "rd_case_reports": 0.04,
        "rd_natural_history": 0.04,
        "rd_registries": 0.02,
        "rd_trials": 0.01,
        "genomic_evidence": 0.01,
    },

    # -- Gene Therapy Eligibility ----------------------------------
    DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY: {
        "rd_therapies": 0.22,
        "rd_trials": 0.15,
        "rd_genes": 0.12,
        "rd_variants": 0.10,
        "rd_guidelines": 0.08,
        "rd_diseases": 0.08,
        "rd_literature": 0.06,
        "rd_case_reports": 0.05,
        "rd_registries": 0.04,
        "rd_natural_history": 0.03,
        "rd_phenotypes": 0.03,
        "rd_pathways": 0.02,
        "genomic_evidence": 0.01,
        "rd_newborn_screening": 0.01,
    },

    # -- Undiagnosed Disease Program -------------------------------
    DiagnosticWorkflowType.UNDIAGNOSED_DISEASE: {
        "rd_phenotypes": 0.15,
        "rd_genes": 0.14,
        "rd_variants": 0.14,
        "rd_diseases": 0.12,
        "rd_case_reports": 0.10,
        "rd_literature": 0.08,
        "rd_pathways": 0.06,
        "rd_guidelines": 0.05,
        "rd_natural_history": 0.04,
        "genomic_evidence": 0.04,
        "rd_therapies": 0.03,
        "rd_trials": 0.02,
        "rd_registries": 0.02,
        "rd_newborn_screening": 0.01,
    },

    # -- General (no specific workflow) ----------------------------
    DiagnosticWorkflowType.GENERAL: {
        "rd_phenotypes": 0.12,
        "rd_diseases": 0.11,
        "rd_genes": 0.10,
        "rd_variants": 0.10,
        "rd_literature": 0.08,
        "rd_case_reports": 0.07,
        "rd_therapies": 0.07,
        "rd_trials": 0.06,
        "rd_guidelines": 0.06,
        "rd_pathways": 0.06,
        "rd_natural_history": 0.05,
        "rd_newborn_screening": 0.05,
        "rd_registries": 0.04,
        "genomic_evidence": 0.03,
    },
}
"""Per-workflow boosted search weights.

Each workflow maps every collection to a weight that sums to ~1.0.
The dominant collection for the workflow receives the highest weight
so that domain-relevant evidence is surfaced preferentially during
multi-collection RAG retrieval.
"""


# ===================================================================
# HELPER FUNCTIONS
# ===================================================================


def get_collection_config(name: str) -> CollectionConfig:
    """Look up a :class:`CollectionConfig` by full collection name.

    Args:
        name: Full Milvus collection name (e.g. ``rd_phenotypes``)
            **or** a short alias (e.g. ``phenotypes``).

    Returns:
        The matching :class:`CollectionConfig`.

    Raises:
        ValueError: If *name* does not match any known collection.
    """
    # Direct lookup by full name
    for cfg in ALL_COLLECTIONS:
        if cfg.name == name:
            return cfg

    # Fallback: resolve short alias first
    resolved = COLLECTION_NAMES.get(name)
    if resolved is not None:
        for cfg in ALL_COLLECTIONS:
            if cfg.name == resolved:
                return cfg

    valid = [cfg.name for cfg in ALL_COLLECTIONS]
    raise ValueError(
        f"Unknown collection '{name}'. "
        f"Valid collections: {valid}"
    )


def get_all_collection_names() -> List[str]:
    """Return a list of all 14 full Milvus collection names.

    Returns:
        Ordered list of collection name strings.
    """
    return [cfg.name for cfg in ALL_COLLECTIONS]


def get_search_weights(
    workflow: Optional[DiagnosticWorkflowType] = None,
) -> Dict[str, float]:
    """Return collection search weights, optionally boosted for a workflow.

    When *workflow* is ``None`` (or not found in the boost table), the
    default base weights from each :class:`CollectionConfig` are returned.

    Args:
        workflow: Optional :class:`DiagnosticWorkflowType` to apply
            workflow-specific weight boosting.

    Returns:
        Dict mapping collection name to its search weight (0.0 - 1.0).
    """
    if workflow is not None and workflow in WORKFLOW_COLLECTION_WEIGHTS:
        return dict(WORKFLOW_COLLECTION_WEIGHTS[workflow])
    return dict(_DEFAULT_SEARCH_WEIGHTS)
