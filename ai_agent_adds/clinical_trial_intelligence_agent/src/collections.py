"""Milvus collection schemas for Clinical Trial Intelligence Agent.

Defines 14 domain-specific vector collections for clinical trial intelligence:
  - trial_protocols          — Clinical trial protocol documents and metadata
  - trial_eligibility        — Inclusion/exclusion eligibility criteria
  - trial_endpoints          — Primary, secondary, and exploratory endpoints
  - trial_sites              — Investigational site data and performance
  - trial_investigators      — Principal investigators and research profiles
  - trial_results            — Published trial results and outcomes
  - trial_regulatory         — Regulatory submissions, decisions, and documents
  - trial_literature         — Published clinical research literature
  - trial_biomarkers         — Biomarker assays, thresholds, and validation
  - trial_safety             — Adverse events and safety signal data
  - trial_rwe               — Real-world evidence studies and registries
  - trial_adaptive           — Adaptive trial designs and decision rules
  - trial_guidelines         — ICH/FDA/EMA clinical trial guidelines
  - genomic_evidence         — Shared genomic evidence (read-only)

Follows the same pymilvus pattern as:
  cardiology_intelligence_agent/src/collections.py

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

from src.models import TrialWorkflowType


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

EMBEDDING_DIM = 384       # BGE-small-en-v1.5
INDEX_TYPE = "IVF_FLAT"
METRIC_TYPE = "COSINE"
NLIST = 128


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION CONFIG DATACLASS
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class CollectionConfig:
    """Configuration for a single Milvus vector collection.

    Attributes:
        name: Milvus collection name (e.g. ``trial_protocols``).
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


# ═══════════════════════════════════════════════════════════════════════
# HELPER — EMBEDDING FIELD
# ═══════════════════════════════════════════════════════════════════════


def _make_embedding_field() -> FieldSchema:
    """Create the standard 384-dim FLOAT_VECTOR embedding field.

    All 14 clinical trial collections share the same embedding specification
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


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION SCHEMA DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

# ── trial_protocols ──────────────────────────────────────────────────

PROTOCOLS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="trial_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Official trial title",
    ),
    FieldSchema(
        name="phase",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Trial phase (Phase I, Phase II, Phase III, Phase IV)",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Trial recruitment status",
    ),
    FieldSchema(
        name="therapeutic_area",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Therapeutic area (oncology, cardiology, neurology, etc.)",
    ),
    FieldSchema(
        name="sponsor",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Lead sponsor organization",
    ),
    FieldSchema(
        name="start_date",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Study start date (YYYY-MM-DD or YYYY-MM)",
    ),
    FieldSchema(
        name="enrollment_target",
        dtype=DataType.INT32,
        description="Target participant enrollment count",
    ),
    FieldSchema(
        name="text_content",
        dtype=DataType.VARCHAR,
        max_length=8192,
        description="Protocol summary or full text content chunk",
    ),
]

PROTOCOLS_CONFIG = CollectionConfig(
    name="trial_protocols",
    description="Clinical trial protocol documents with phase, sponsor, and enrollment metadata",
    schema_fields=PROTOCOLS_FIELDS,
    estimated_records=5000,
    search_weight=0.10,
)

# ── trial_eligibility ────────────────────────────────────────────────

ELIGIBILITY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="trial_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="criterion_type",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Criterion type (inclusion or exclusion)",
    ),
    FieldSchema(
        name="criterion_text",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Full eligibility criterion text",
    ),
    FieldSchema(
        name="logic_operator",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Logic operator (AND, OR, NOT) for compound criteria",
    ),
    FieldSchema(
        name="population_impact",
        dtype=DataType.FLOAT,
        description="Estimated population impact (fraction excluded, 0.0 - 1.0)",
    ),
]

ELIGIBILITY_CONFIG = CollectionConfig(
    name="trial_eligibility",
    description="Trial inclusion and exclusion eligibility criteria with population impact estimates",
    schema_fields=ELIGIBILITY_FIELDS,
    estimated_records=50000,
    search_weight=0.09,
)

# ── trial_endpoints ──────────────────────────────────────────────────

ENDPOINTS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="trial_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="endpoint_type",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Endpoint type (primary, secondary, exploratory, safety)",
    ),
    FieldSchema(
        name="measure",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Outcome measure description",
    ),
    FieldSchema(
        name="time_frame",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Time frame for endpoint assessment",
    ),
    FieldSchema(
        name="statistical_method",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Statistical analysis method (e.g., log-rank, MMRM, Cox PH)",
    ),
]

ENDPOINTS_CONFIG = CollectionConfig(
    name="trial_endpoints",
    description="Clinical trial endpoint definitions with measures, time frames, and statistical methods",
    schema_fields=ENDPOINTS_FIELDS,
    estimated_records=20000,
    search_weight=0.08,
)

# ── trial_sites ──────────────────────────────────────────────────────

SITES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="trial_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="site_id",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Site identifier",
    ),
    FieldSchema(
        name="facility_name",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Facility or institution name",
    ),
    FieldSchema(
        name="city",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="City",
    ),
    FieldSchema(
        name="state",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="State or province",
    ),
    FieldSchema(
        name="country",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Country",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Site recruitment status",
    ),
    FieldSchema(
        name="enrollment_count",
        dtype=DataType.INT32,
        description="Number of participants enrolled at this site",
    ),
]

SITES_CONFIG = CollectionConfig(
    name="trial_sites",
    description="Investigational site locations, enrollment counts, and recruitment status",
    schema_fields=SITES_FIELDS,
    estimated_records=30000,
    search_weight=0.07,
)

# ── trial_investigators ──────────────────────────────────────────────

INVESTIGATORS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="investigator_id",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Investigator unique identifier",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Investigator full name",
    ),
    FieldSchema(
        name="specialty",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Medical specialty",
    ),
    FieldSchema(
        name="h_index",
        dtype=DataType.INT32,
        description="H-index publication metric",
    ),
    FieldSchema(
        name="publication_count",
        dtype=DataType.INT32,
        description="Total publication count",
    ),
    FieldSchema(
        name="therapeutic_areas",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Semicolon-separated therapeutic area expertise",
    ),
]

INVESTIGATORS_CONFIG = CollectionConfig(
    name="trial_investigators",
    description="Principal investigators with publication metrics and therapeutic area expertise",
    schema_fields=INVESTIGATORS_FIELDS,
    estimated_records=5000,
    search_weight=0.05,
)

# ── trial_results ────────────────────────────────────────────────────

RESULTS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="trial_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="outcome",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Outcome description and result summary",
    ),
    FieldSchema(
        name="p_value",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Statistical p-value",
    ),
    FieldSchema(
        name="effect_size",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Effect size (HR, OR, RR, mean difference, etc.)",
    ),
    FieldSchema(
        name="confidence_interval",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="95% confidence interval",
    ),
    FieldSchema(
        name="publication_pmid",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="PubMed identifier of the results publication",
    ),
]

RESULTS_CONFIG = CollectionConfig(
    name="trial_results",
    description="Published trial results with statistical outcomes and effect sizes",
    schema_fields=RESULTS_FIELDS,
    estimated_records=3000,
    search_weight=0.09,
)

# ── trial_regulatory ─────────────────────────────────────────────────

REGULATORY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="submission_id",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Regulatory submission identifier",
    ),
    FieldSchema(
        name="agency",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Regulatory agency (FDA, EMA, PMDA, etc.)",
    ),
    FieldSchema(
        name="decision",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Regulatory decision (approved, refused, pending, withdrawn)",
    ),
    FieldSchema(
        name="document_type",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Document type (IND, CSR, briefing, PSP, RMP, DSUR)",
    ),
    FieldSchema(
        name="drug_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Drug or product name",
    ),
    FieldSchema(
        name="indication",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Approved or proposed indication",
    ),
]

REGULATORY_CONFIG = CollectionConfig(
    name="trial_regulatory",
    description="Regulatory submissions, decisions, and documentation across agencies",
    schema_fields=REGULATORY_FIELDS,
    estimated_records=2000,
    search_weight=0.07,
)

# ── trial_literature ─────────────────────────────────────────────────

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
        max_length=32,
        description="PubMed identifier",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Publication title",
    ),
    FieldSchema(
        name="journal",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Journal name",
    ),
    FieldSchema(
        name="mesh_terms",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Semicolon-separated MeSH terms",
    ),
    FieldSchema(
        name="publication_year",
        dtype=DataType.INT32,
        description="Publication year",
    ),
    FieldSchema(
        name="study_type",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Study design (RCT, meta-analysis, cohort, systematic review)",
    ),
]

LITERATURE_CONFIG = CollectionConfig(
    name="trial_literature",
    description="Published clinical trial methodology and results literature",
    schema_fields=LITERATURE_FIELDS,
    estimated_records=10000,
    search_weight=0.08,
)

# ── trial_biomarkers ─────────────────────────────────────────────────

BIOMARKERS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="biomarker",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Biomarker name (e.g., PD-L1, HER2, BRCA1, TMB)",
    ),
    FieldSchema(
        name="assay",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Assay or testing platform",
    ),
    FieldSchema(
        name="threshold",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Positivity or actionability threshold",
    ),
    FieldSchema(
        name="validated",
        dtype=DataType.BOOL,
        description="Whether the biomarker is analytically validated",
    ),
    FieldSchema(
        name="trial_context",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Trial context for biomarker use (stratification, selection, endpoint)",
    ),
]

BIOMARKERS_CONFIG = CollectionConfig(
    name="trial_biomarkers",
    description="Biomarker assays, thresholds, and validation status in trial contexts",
    schema_fields=BIOMARKERS_FIELDS,
    estimated_records=3000,
    search_weight=0.07,
)

# ── trial_safety ─────────────────────────────────────────────────────

SAFETY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="trial_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="event_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Adverse event MedDRA preferred term",
    ),
    FieldSchema(
        name="severity",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="CTCAE severity grade or classification",
    ),
    FieldSchema(
        name="frequency",
        dtype=DataType.FLOAT,
        description="Observed frequency in the study population (0.0 - 1.0)",
    ),
    FieldSchema(
        name="soc_term",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="MedDRA System Organ Class term",
    ),
]

SAFETY_CONFIG = CollectionConfig(
    name="trial_safety",
    description="Adverse event profiles with MedDRA coding, severity, and frequency data",
    schema_fields=SAFETY_FIELDS,
    estimated_records=20000,
    search_weight=0.08,
)

# ── trial_rwe ────────────────────────────────────────────────────────

RWE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="RWE data source (e.g., claims database, EHR registry, patient registry)",
    ),
    FieldSchema(
        name="population",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Study population description",
    ),
    FieldSchema(
        name="outcome",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Real-world outcome or effectiveness measure",
    ),
    FieldSchema(
        name="study_design",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Study design (retrospective cohort, prospective registry, etc.)",
    ),
    FieldSchema(
        name="sample_size",
        dtype=DataType.INT32,
        description="Study sample size",
    ),
]

RWE_CONFIG = CollectionConfig(
    name="trial_rwe",
    description="Real-world evidence studies from claims, EHR, and patient registries",
    schema_fields=RWE_FIELDS,
    estimated_records=2000,
    search_weight=0.06,
)

# ── trial_adaptive ───────────────────────────────────────────────────

ADAPTIVE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="design_type",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Adaptive design type (Bayesian, seamless, platform, basket, umbrella)",
    ),
    FieldSchema(
        name="decision_rule",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Interim decision rule (e.g., futility, efficacy, dose modification)",
    ),
    FieldSchema(
        name="trigger_criteria",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Criteria that trigger adaptation (sample size, information fraction)",
    ),
    FieldSchema(
        name="historical_precedent",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Historical precedent and regulatory acceptance of this design",
    ),
]

ADAPTIVE_CONFIG = CollectionConfig(
    name="trial_adaptive",
    description="Adaptive trial designs, decision rules, and regulatory precedents",
    schema_fields=ADAPTIVE_FIELDS,
    estimated_records=500,
    search_weight=0.05,
)

# ── trial_guidelines ─────────────────────────────────────────────────

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
        description="Guideline identifier (e.g., ICH E6(R3), ICH E9(R1))",
    ),
    FieldSchema(
        name="organization",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Issuing organization (ICH, FDA, EMA, WHO)",
    ),
    FieldSchema(
        name="version",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Guideline version or revision",
    ),
    FieldSchema(
        name="recommendation_text",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Guideline recommendation or requirement text",
    ),
    FieldSchema(
        name="evidence_class",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Evidence class supporting the recommendation",
    ),
]

GUIDELINES_CONFIG = CollectionConfig(
    name="trial_guidelines",
    description="ICH, FDA, and EMA clinical trial guidelines and regulatory requirements",
    schema_fields=GUIDELINES_FIELDS,
    estimated_records=1000,
    search_weight=0.08,
)

# ── genomic_evidence (shared, read-only) ─────────────────────────────

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
        max_length=64,
        description="Gene symbol (HUGO nomenclature)",
    ),
    FieldSchema(
        name="variant",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Variant designation (e.g. rs ID, HGVS notation)",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="ClinVar clinical significance (pathogenic, likely pathogenic, VUS, etc.)",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Associated condition or disease",
    ),
    FieldSchema(
        name="evidence_summary",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Evidence summary text",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Data source (ClinVar, AlphaMissense, gnomAD)",
    ),
]

GENOMIC_CONFIG = CollectionConfig(
    name="genomic_evidence",
    description="Shared genomic variant evidence from ClinVar, AlphaMissense, and gnomAD (read-only)",
    schema_fields=GENOMIC_FIELDS,
    estimated_records=100000,
    search_weight=0.03,
)


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION REGISTRY
# ═══════════════════════════════════════════════════════════════════════

ALL_COLLECTIONS: List[CollectionConfig] = [
    PROTOCOLS_CONFIG,
    ELIGIBILITY_CONFIG,
    ENDPOINTS_CONFIG,
    SITES_CONFIG,
    INVESTIGATORS_CONFIG,
    RESULTS_CONFIG,
    REGULATORY_CONFIG,
    LITERATURE_CONFIG,
    BIOMARKERS_CONFIG,
    SAFETY_CONFIG,
    RWE_CONFIG,
    ADAPTIVE_CONFIG,
    GUIDELINES_CONFIG,
    GENOMIC_CONFIG,
]
"""Ordered list of all 14 clinical trial collection configurations."""


COLLECTION_NAMES: Dict[str, str] = {
    "protocols": "trial_protocols",
    "eligibility": "trial_eligibility",
    "endpoints": "trial_endpoints",
    "sites": "trial_sites",
    "investigators": "trial_investigators",
    "results": "trial_results",
    "regulatory": "trial_regulatory",
    "literature": "trial_literature",
    "biomarkers": "trial_biomarkers",
    "safety": "trial_safety",
    "rwe": "trial_rwe",
    "adaptive": "trial_adaptive",
    "guidelines": "trial_guidelines",
    "genomic": "genomic_evidence",
}
"""Mapping of short alias names to full Milvus collection names."""


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION SCHEMAS (pymilvus CollectionSchema objects)
# ═══════════════════════════════════════════════════════════════════════

COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    cfg.name: CollectionSchema(
        fields=cfg.schema_fields,
        description=cfg.description,
    )
    for cfg in ALL_COLLECTIONS
}
"""Mapping of collection name to pymilvus CollectionSchema, ready for
``Collection(name=..., schema=...)`` instantiation."""


# ═══════════════════════════════════════════════════════════════════════
# DEFAULT SEARCH WEIGHTS
# ═══════════════════════════════════════════════════════════════════════

_DEFAULT_SEARCH_WEIGHTS: Dict[str, float] = {
    cfg.name: cfg.search_weight for cfg in ALL_COLLECTIONS
}
"""Base search weights used when no workflow-specific boost is applied.
Sum: {sum:.2f}.""".format(sum=sum(cfg.search_weight for cfg in ALL_COLLECTIONS))


# ═══════════════════════════════════════════════════════════════════════
# WORKFLOW-SPECIFIC COLLECTION WEIGHTS
# ═══════════════════════════════════════════════════════════════════════

WORKFLOW_COLLECTION_WEIGHTS: Dict[TrialWorkflowType, Dict[str, float]] = {
    # ── Protocol Design ──────────────────────────────────────────────
    TrialWorkflowType.PROTOCOL_DESIGN: {
        "trial_protocols": 0.20,
        "trial_endpoints": 0.15,
        "trial_eligibility": 0.12,
        "trial_guidelines": 0.10,
        "trial_adaptive": 0.08,
        "trial_literature": 0.08,
        "trial_results": 0.06,
        "trial_regulatory": 0.05,
        "trial_biomarkers": 0.04,
        "trial_safety": 0.04,
        "trial_sites": 0.03,
        "trial_rwe": 0.02,
        "trial_investigators": 0.02,
        "genomic_evidence": 0.01,
    },

    # ── Patient Matching ─────────────────────────────────────────────
    TrialWorkflowType.PATIENT_MATCHING: {
        "trial_eligibility": 0.25,
        "trial_protocols": 0.15,
        "trial_sites": 0.12,
        "trial_biomarkers": 0.10,
        "genomic_evidence": 0.08,
        "trial_endpoints": 0.06,
        "trial_safety": 0.05,
        "trial_results": 0.05,
        "trial_guidelines": 0.04,
        "trial_literature": 0.03,
        "trial_investigators": 0.03,
        "trial_rwe": 0.02,
        "trial_regulatory": 0.01,
        "trial_adaptive": 0.01,
    },

    # ── Site Selection ───────────────────────────────────────────────
    TrialWorkflowType.SITE_SELECTION: {
        "trial_sites": 0.25,
        "trial_investigators": 0.18,
        "trial_protocols": 0.10,
        "trial_results": 0.08,
        "trial_eligibility": 0.08,
        "trial_rwe": 0.06,
        "trial_guidelines": 0.05,
        "trial_literature": 0.05,
        "trial_regulatory": 0.04,
        "trial_safety": 0.04,
        "trial_endpoints": 0.03,
        "trial_biomarkers": 0.02,
        "trial_adaptive": 0.01,
        "genomic_evidence": 0.01,
    },

    # ── Eligibility Optimization ─────────────────────────────────────
    TrialWorkflowType.ELIGIBILITY_OPTIMIZATION: {
        "trial_eligibility": 0.25,
        "trial_protocols": 0.12,
        "trial_rwe": 0.10,
        "trial_results": 0.10,
        "trial_guidelines": 0.08,
        "trial_literature": 0.08,
        "trial_biomarkers": 0.06,
        "trial_safety": 0.05,
        "trial_sites": 0.04,
        "trial_endpoints": 0.04,
        "genomic_evidence": 0.03,
        "trial_regulatory": 0.02,
        "trial_investigators": 0.02,
        "trial_adaptive": 0.01,
    },

    # ── Adaptive Design ──────────────────────────────────────────────
    TrialWorkflowType.ADAPTIVE_DESIGN: {
        "trial_adaptive": 0.25,
        "trial_endpoints": 0.15,
        "trial_guidelines": 0.12,
        "trial_protocols": 0.10,
        "trial_regulatory": 0.08,
        "trial_results": 0.06,
        "trial_literature": 0.06,
        "trial_biomarkers": 0.05,
        "trial_safety": 0.04,
        "trial_eligibility": 0.03,
        "trial_sites": 0.02,
        "trial_rwe": 0.02,
        "trial_investigators": 0.01,
        "genomic_evidence": 0.01,
    },

    # ── Safety Signal Detection ──────────────────────────────────────
    TrialWorkflowType.SAFETY_SIGNAL: {
        "trial_safety": 0.25,
        "trial_results": 0.15,
        "trial_protocols": 0.10,
        "trial_literature": 0.10,
        "trial_rwe": 0.08,
        "trial_guidelines": 0.06,
        "trial_regulatory": 0.06,
        "trial_biomarkers": 0.05,
        "trial_endpoints": 0.04,
        "trial_eligibility": 0.04,
        "trial_adaptive": 0.03,
        "trial_sites": 0.02,
        "trial_investigators": 0.01,
        "genomic_evidence": 0.01,
    },

    # ── Regulatory Documents ─────────────────────────────────────────
    TrialWorkflowType.REGULATORY_DOCS: {
        "trial_regulatory": 0.25,
        "trial_guidelines": 0.18,
        "trial_results": 0.12,
        "trial_safety": 0.10,
        "trial_protocols": 0.08,
        "trial_literature": 0.06,
        "trial_endpoints": 0.05,
        "trial_biomarkers": 0.04,
        "trial_adaptive": 0.03,
        "trial_eligibility": 0.03,
        "trial_rwe": 0.03,
        "trial_sites": 0.01,
        "trial_investigators": 0.01,
        "genomic_evidence": 0.01,
    },

    # ── Competitive Intelligence ─────────────────────────────────────
    TrialWorkflowType.COMPETITIVE_INTEL: {
        "trial_protocols": 0.20,
        "trial_results": 0.15,
        "trial_endpoints": 0.12,
        "trial_eligibility": 0.10,
        "trial_literature": 0.08,
        "trial_sites": 0.08,
        "trial_regulatory": 0.06,
        "trial_biomarkers": 0.05,
        "trial_safety": 0.05,
        "trial_investigators": 0.04,
        "trial_adaptive": 0.03,
        "trial_rwe": 0.02,
        "trial_guidelines": 0.01,
        "genomic_evidence": 0.01,
    },

    # ── Diversity Assessment ─────────────────────────────────────────
    TrialWorkflowType.DIVERSITY_ASSESSMENT: {
        "trial_sites": 0.22,
        "trial_eligibility": 0.18,
        "trial_protocols": 0.12,
        "trial_rwe": 0.10,
        "trial_guidelines": 0.08,
        "trial_results": 0.06,
        "trial_literature": 0.06,
        "trial_investigators": 0.05,
        "trial_biomarkers": 0.04,
        "trial_safety": 0.03,
        "trial_endpoints": 0.02,
        "trial_regulatory": 0.02,
        "trial_adaptive": 0.01,
        "genomic_evidence": 0.01,
    },

    # ── Decentralized Planning ───────────────────────────────────────
    TrialWorkflowType.DECENTRALIZED_PLANNING: {
        "trial_sites": 0.18,
        "trial_protocols": 0.15,
        "trial_guidelines": 0.12,
        "trial_endpoints": 0.10,
        "trial_eligibility": 0.08,
        "trial_regulatory": 0.08,
        "trial_literature": 0.06,
        "trial_biomarkers": 0.05,
        "trial_safety": 0.05,
        "trial_results": 0.04,
        "trial_rwe": 0.04,
        "trial_adaptive": 0.02,
        "trial_investigators": 0.02,
        "genomic_evidence": 0.01,
    },

    # ── General (no specific workflow) ───────────────────────────────
    TrialWorkflowType.GENERAL: {
        "trial_protocols": 0.10,
        "trial_eligibility": 0.09,
        "trial_endpoints": 0.08,
        "trial_results": 0.09,
        "trial_safety": 0.08,
        "trial_literature": 0.08,
        "trial_guidelines": 0.08,
        "trial_sites": 0.07,
        "trial_regulatory": 0.07,
        "trial_biomarkers": 0.07,
        "trial_rwe": 0.06,
        "trial_adaptive": 0.05,
        "trial_investigators": 0.05,
        "genomic_evidence": 0.03,
    },
}
"""Per-workflow boosted search weights.

Each workflow maps every collection to a weight that sums to ~1.0.
The dominant collection for the workflow receives the highest weight
so that domain-relevant evidence is surfaced preferentially during
multi-collection RAG retrieval.
"""


# ═══════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════


def get_collection_config(name: str) -> CollectionConfig:
    """Look up a :class:`CollectionConfig` by full collection name.

    Args:
        name: Full Milvus collection name (e.g. ``trial_protocols``)
            **or** a short alias (e.g. ``protocols``).

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
    workflow: Optional[TrialWorkflowType] = None,
) -> Dict[str, float]:
    """Return collection search weights, optionally boosted for a workflow.

    When *workflow* is ``None`` (or not found in the boost table), the
    default base weights from each :class:`CollectionConfig` are returned.

    Args:
        workflow: Optional :class:`TrialWorkflowType` to apply
            workflow-specific weight boosting.

    Returns:
        Dict mapping collection name to its search weight (0.0 - 1.0).
    """
    if workflow is not None and workflow in WORKFLOW_COLLECTION_WEIGHTS:
        return dict(WORKFLOW_COLLECTION_WEIGHTS[workflow])
    return dict(_DEFAULT_SEARCH_WEIGHTS)
