"""Milvus collection management for Biomarker Intelligence Agent.

Manages 11 biomarker collections (10 domain-specific + 1 read-only genomic):
  - biomarker_reference          — Reference biomarker definitions & ranges
  - biomarker_genetic_variants   — Genetic variants affecting biomarkers
  - biomarker_pgx_rules          — Pharmacogenomic dosing rules (CPIC)
  - biomarker_disease_trajectories — Disease progression trajectories
  - biomarker_clinical_evidence  — Published clinical evidence
  - biomarker_nutrition          — Genotype-aware nutrition guidelines
  - biomarker_drug_interactions  — Gene-drug interactions
  - biomarker_aging_markers      — Epigenetic aging clock markers
  - biomarker_genotype_adjustments — Genotype-based reference range adjustments
  - biomarker_monitoring         — Condition-specific monitoring protocols
  + genomic_evidence             — Shared genomic variants (read-only)

Follows the same pymilvus pattern as:
  cart_intelligence_agent/src/collections.py (CARTCollectionManager)
  rag-chat-pipeline/src/milvus_client.py (MilvusClient)

Author: Adam Jones
Date: March 2026
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

from loguru import logger
from pymilvus import (
    Collection,
    CollectionSchema,
    DataType,
    FieldSchema,
    connections,
    utility,
)

from src.models import (
    AgingMarker,
    AJCarrierScreeningEntry,
    BiomarkerReference,
    ClinicalEvidence,
    CriticalValue,
    DiscordanceRule,
    DiseaseTrajectory,
    DrugInteraction,
    GeneticVariant,
    GenotypeAdjustment,
    MonitoringProtocol,
    NutritionGuideline,
    PGxRule,
)

# ═══════════════════════════════════════════════════════════════════════
# COLLECTION SCHEMA DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

EMBEDDING_DIM = 384  # BGE-small-en-v1.5

# ── biomarker_reference ─────────────────────────────────────────────

BIOMARKER_REFERENCE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique biomarker identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Biomarker display name",
    ),
    FieldSchema(
        name="unit",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Measurement unit (e.g., mg/dL)",
    ),
    FieldSchema(
        name="category",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="CBC, CMP, Lipids, Thyroid, Inflammation, Nutrients",
    ),
    FieldSchema(
        name="ref_range_min",
        dtype=DataType.FLOAT,
        description="Standard reference range lower bound",
    ),
    FieldSchema(
        name="ref_range_max",
        dtype=DataType.FLOAT,
        description="Standard reference range upper bound",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Clinical interpretation and significance",
    ),
    FieldSchema(
        name="epigenetic_clock",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="PhenoAge or GrimAge coefficient if applicable",
    ),
    FieldSchema(
        name="genetic_modifiers",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated genes that modify this biomarker",
    ),
]

BIOMARKER_REFERENCE_SCHEMA = CollectionSchema(
    fields=BIOMARKER_REFERENCE_FIELDS,
    description="Reference biomarker definitions, ranges, and clinical significance",
)

# ── biomarker_genetic_variants ──────────────────────────────────────

GENETIC_VARIANTS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique variant identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Gene symbol (e.g., MTHFR)",
    ),
    FieldSchema(
        name="rs_id",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="dbSNP rsID (e.g., rs1801133)",
    ),
    FieldSchema(
        name="risk_allele",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Risk allele",
    ),
    FieldSchema(
        name="protective_allele",
        dtype=DataType.VARCHAR,
        max_length=5,
        description="Protective allele",
    ),
    FieldSchema(
        name="effect_size",
        dtype=DataType.VARCHAR,
        max_length=250,
        description="Effect size description",
    ),
    FieldSchema(
        name="mechanism",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Molecular mechanism of the variant effect",
    ),
    FieldSchema(
        name="disease_associations",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Comma-separated disease associations",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

GENETIC_VARIANTS_SCHEMA = CollectionSchema(
    fields=GENETIC_VARIANTS_FIELDS,
    description="Genetic variants affecting biomarker levels and disease risk",
)

# ── biomarker_pgx_rules ────────────────────────────────────────────

PGX_RULES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique PGx rule identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Pharmacogene (e.g., CYP2D6)",
    ),
    FieldSchema(
        name="star_alleles",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Star allele combination (e.g., *1/*2)",
    ),
    FieldSchema(
        name="drug",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Drug name",
    ),
    FieldSchema(
        name="phenotype",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Metabolizer phenotype (ultra_rapid, normal, intermediate, poor)",
    ),
    FieldSchema(
        name="cpic_level",
        dtype=DataType.VARCHAR,
        max_length=10,
        description="CPIC evidence level (1A, 1B, 2A, 2B, 3)",
    ),
    FieldSchema(
        name="recommendation",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="CPIC dosing recommendation",
    ),
    FieldSchema(
        name="evidence_url",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="URL to CPIC guideline or PharmGKB entry",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

PGX_RULES_SCHEMA = CollectionSchema(
    fields=PGX_RULES_FIELDS,
    description="Pharmacogenomic dosing rules from CPIC guidelines",
)

# ── biomarker_disease_trajectories ──────────────────────────────────

DISEASE_TRAJECTORIES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique trajectory identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Disease category",
    ),
    FieldSchema(
        name="stage",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Disease stage (e.g., pre-diabetes, early, advanced)",
    ),
    FieldSchema(
        name="biomarker_pattern",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="JSON string of biomarker thresholds defining this stage",
    ),
    FieldSchema(
        name="years_to_diagnosis",
        dtype=DataType.FLOAT,
        description="Estimated years from this stage to clinical diagnosis",
    ),
    FieldSchema(
        name="intervention_window",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Description of the intervention opportunity",
    ),
    FieldSchema(
        name="risk_reduction_pct",
        dtype=DataType.FLOAT,
        description="Potential risk reduction percentage with intervention",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

DISEASE_TRAJECTORIES_SCHEMA = CollectionSchema(
    fields=DISEASE_TRAJECTORIES_FIELDS,
    description="Disease progression trajectories with biomarker patterns",
)

# ── biomarker_clinical_evidence ─────────────────────────────────────

CLINICAL_EVIDENCE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique evidence identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="pmid",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="PubMed ID",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Publication title",
    ),
    FieldSchema(
        name="finding",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Key finding from the publication",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Publication year",
    ),
    FieldSchema(
        name="disease_area",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Disease area or specialty",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

CLINICAL_EVIDENCE_SCHEMA = CollectionSchema(
    fields=CLINICAL_EVIDENCE_FIELDS,
    description="Published clinical evidence for biomarker interpretation",
)

# ── biomarker_nutrition ─────────────────────────────────────────────

NUTRITION_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique guideline identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="nutrient",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Nutrient name (e.g., Folate)",
    ),
    FieldSchema(
        name="genetic_context",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Genetic context (e.g., MTHFR C677T heterozygous)",
    ),
    FieldSchema(
        name="recommended_form",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Recommended supplement form (e.g., methylfolate)",
    ),
    FieldSchema(
        name="dose_range",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Dosing range (e.g., 400-800 mcg/day)",
    ),
    FieldSchema(
        name="evidence_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Summary of evidence supporting this guideline",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

NUTRITION_SCHEMA = CollectionSchema(
    fields=NUTRITION_FIELDS,
    description="Genotype-aware nutrition guidelines and supplement recommendations",
)

# ── biomarker_drug_interactions ─────────────────────────────────────

DRUG_INTERACTIONS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique interaction identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="drug",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Drug name",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Gene involved (e.g., CYP2C19)",
    ),
    FieldSchema(
        name="interaction_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Interaction type (substrate, inhibitor, inducer)",
    ),
    FieldSchema(
        name="severity",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Severity (major, moderate, minor)",
    ),
    FieldSchema(
        name="alternative",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Alternative drug recommendation",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

DRUG_INTERACTIONS_SCHEMA = CollectionSchema(
    fields=DRUG_INTERACTIONS_FIELDS,
    description="Gene-drug interactions for pharmacogenomic-aware prescribing",
)

# ── biomarker_aging_markers ─────────────────────────────────────────

AGING_MARKERS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique marker identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="marker_name",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Marker name (e.g., Albumin)",
    ),
    FieldSchema(
        name="clock_type",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="PhenoAge or GrimAge",
    ),
    FieldSchema(
        name="coefficient",
        dtype=DataType.FLOAT,
        description="Coefficient weight in the aging clock algorithm",
    ),
    FieldSchema(
        name="unit",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Measurement unit",
    ),
    FieldSchema(
        name="interpretation",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Clinical interpretation in aging context",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

AGING_MARKERS_SCHEMA = CollectionSchema(
    fields=AGING_MARKERS_FIELDS,
    description="Epigenetic aging clock markers (PhenoAge, GrimAge)",
)

# ── biomarker_genotype_adjustments ──────────────────────────────────

GENOTYPE_ADJUSTMENTS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique adjustment identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="biomarker",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Biomarker name",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Gene symbol",
    ),
    FieldSchema(
        name="rs_id",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="dbSNP rsID",
    ),
    FieldSchema(
        name="genotype_ref",
        dtype=DataType.VARCHAR,
        max_length=10,
        description="Reference genotype (e.g., CC)",
    ),
    FieldSchema(
        name="genotype_het",
        dtype=DataType.VARCHAR,
        max_length=10,
        description="Heterozygous genotype (e.g., CT)",
    ),
    FieldSchema(
        name="genotype_hom",
        dtype=DataType.VARCHAR,
        max_length=10,
        description="Homozygous alternate genotype (e.g., TT)",
    ),
    FieldSchema(
        name="adjusted_min",
        dtype=DataType.FLOAT,
        description="Adjusted reference range lower bound",
    ),
    FieldSchema(
        name="adjusted_max",
        dtype=DataType.FLOAT,
        description="Adjusted reference range upper bound",
    ),
    FieldSchema(
        name="rationale",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Rationale for the genotype-based adjustment",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

GENOTYPE_ADJUSTMENTS_SCHEMA = CollectionSchema(
    fields=GENOTYPE_ADJUSTMENTS_FIELDS,
    description="Genotype-based reference range adjustments for biomarkers",
)

# ── biomarker_monitoring ────────────────────────────────────────────

MONITORING_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique protocol identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Medical condition",
    ),
    FieldSchema(
        name="biomarker_panel",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated biomarker panel",
    ),
    FieldSchema(
        name="frequency",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Monitoring frequency (e.g., every 3 months)",
    ),
    FieldSchema(
        name="trigger_values",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Threshold values that trigger clinical action",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

MONITORING_SCHEMA = CollectionSchema(
    fields=MONITORING_FIELDS,
    description="Condition-specific biomarker monitoring protocols",
)

# ── biomarker_critical_values ──────────────────────────────────────

CRITICAL_VALUES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique critical value identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="biomarker",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Biomarker name",
    ),
    FieldSchema(
        name="loinc_code",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="LOINC code",
    ),
    FieldSchema(
        name="critical_low",
        dtype=DataType.FLOAT,
        description="Critical low threshold",
    ),
    FieldSchema(
        name="critical_high",
        dtype=DataType.FLOAT,
        description="Critical high threshold",
    ),
    FieldSchema(
        name="severity",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Severity level (critical/urgent/warning)",
    ),
    FieldSchema(
        name="escalation_target",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Escalation target",
    ),
    FieldSchema(
        name="clinical_action",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Clinical action to take",
    ),
    FieldSchema(
        name="cross_checks",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated cross-check biomarkers",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

CRITICAL_VALUES_SCHEMA = CollectionSchema(
    fields=CRITICAL_VALUES_FIELDS,
    description="Critical biomarker value thresholds and escalation rules",
)

# ── biomarker_discordance_rules ───────────────────────────────────

DISCORDANCE_RULES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique discordance rule identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Discordance rule name",
    ),
    FieldSchema(
        name="biomarker_a",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="First biomarker in discordance pair",
    ),
    FieldSchema(
        name="biomarker_b",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Second biomarker in discordance pair",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Discordance condition pattern",
    ),
    FieldSchema(
        name="expected_relationship",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Expected relationship between the biomarkers",
    ),
    FieldSchema(
        name="differential_diagnosis",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Comma-separated differential diagnoses",
    ),
    FieldSchema(
        name="agent_handoff",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated agent names for handoff",
    ),
    FieldSchema(
        name="priority",
        dtype=DataType.VARCHAR,
        max_length=10,
        description="Priority level (high/medium/low)",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

DISCORDANCE_RULES_SCHEMA = CollectionSchema(
    fields=DISCORDANCE_RULES_FIELDS,
    description="Cross-biomarker discordance patterns and differential diagnoses",
)

# ── biomarker_aj_carrier_screening ────────────────────────────────

AJ_CARRIER_SCREENING_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Unique AJ carrier screening entry identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Gene symbol",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Disease name",
    ),
    FieldSchema(
        name="inheritance",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Inheritance pattern",
    ),
    FieldSchema(
        name="aj_carrier_frequency",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="AJ population carrier frequency",
    ),
    FieldSchema(
        name="general_carrier_frequency",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="General population carrier frequency",
    ),
    FieldSchema(
        name="common_mutations",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated common mutations",
    ),
    FieldSchema(
        name="loinc_code",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="LOINC code",
    ),
    FieldSchema(
        name="method",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Testing method",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Clinical significance",
    ),
    FieldSchema(
        name="reproductive_implications",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Reproductive implications",
    ),
    FieldSchema(
        name="compound_risks",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Comma-separated compound risks",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

AJ_CARRIER_SCREENING_SCHEMA = CollectionSchema(
    fields=AJ_CARRIER_SCREENING_FIELDS,
    description="Ashkenazi Jewish carrier screening genes and mutations",
)

# ── Genomic Evidence (read-only, created by rag-chat-pipeline) ──────

GENOMIC_EVIDENCE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=200),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="chrom", dtype=DataType.VARCHAR, max_length=10),
    FieldSchema(name="pos", dtype=DataType.INT64),
    FieldSchema(name="ref", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="alt", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="qual", dtype=DataType.FLOAT),
    FieldSchema(name="gene", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="consequence", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="impact", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="genotype", dtype=DataType.VARCHAR, max_length=10),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=2000),
    FieldSchema(name="clinical_significance", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="rsid", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="disease_associations", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="am_pathogenicity", dtype=DataType.FLOAT),
    FieldSchema(name="am_class", dtype=DataType.VARCHAR, max_length=30),
]

GENOMIC_EVIDENCE_SCHEMA = CollectionSchema(
    fields=GENOMIC_EVIDENCE_FIELDS,
    description="Genomic variant evidence (read-only, from rag-chat-pipeline)",
)


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION REGISTRY
# ═══════════════════════════════════════════════════════════════════════

COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    "biomarker_reference": BIOMARKER_REFERENCE_SCHEMA,
    "biomarker_genetic_variants": GENETIC_VARIANTS_SCHEMA,
    "biomarker_pgx_rules": PGX_RULES_SCHEMA,
    "biomarker_disease_trajectories": DISEASE_TRAJECTORIES_SCHEMA,
    "biomarker_clinical_evidence": CLINICAL_EVIDENCE_SCHEMA,
    "biomarker_nutrition": NUTRITION_SCHEMA,
    "biomarker_drug_interactions": DRUG_INTERACTIONS_SCHEMA,
    "biomarker_aging_markers": AGING_MARKERS_SCHEMA,
    "biomarker_genotype_adjustments": GENOTYPE_ADJUSTMENTS_SCHEMA,
    "biomarker_monitoring": MONITORING_SCHEMA,
    "biomarker_critical_values": CRITICAL_VALUES_SCHEMA,
    "biomarker_discordance_rules": DISCORDANCE_RULES_SCHEMA,
    "biomarker_aj_carrier_screening": AJ_CARRIER_SCREENING_SCHEMA,
    "genomic_evidence": GENOMIC_EVIDENCE_SCHEMA,
}

# Maps collection names to their Pydantic model class for validation
# genomic_evidence is None because it's read-only (no inserts from this agent)
COLLECTION_MODELS: Dict[str, type] = {
    "biomarker_reference": BiomarkerReference,
    "biomarker_genetic_variants": GeneticVariant,
    "biomarker_pgx_rules": PGxRule,
    "biomarker_disease_trajectories": DiseaseTrajectory,
    "biomarker_clinical_evidence": ClinicalEvidence,
    "biomarker_nutrition": NutritionGuideline,
    "biomarker_drug_interactions": DrugInteraction,
    "biomarker_aging_markers": AgingMarker,
    "biomarker_genotype_adjustments": GenotypeAdjustment,
    "biomarker_monitoring": MonitoringProtocol,
    "biomarker_critical_values": CriticalValue,
    "biomarker_discordance_rules": DiscordanceRule,
    "biomarker_aj_carrier_screening": AJCarrierScreeningEntry,
    "genomic_evidence": None,
}


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MANAGER
# ═══════════════════════════════════════════════════════════════════════


class BiomarkerCollectionManager:
    """Manages 11 Biomarker Milvus collections (10 owned + 1 read-only genomic).

    Provides create/drop/insert/search operations across the full set of
    biomarker domain collections, following the same pymilvus patterns as
    cart_intelligence_agent/src/collections.py (CARTCollectionManager) and
    rag-chat-pipeline/src/milvus_client.py.

    Usage:
        manager = BiomarkerCollectionManager()
        manager.connect()
        manager.create_all_collections()
        stats = manager.get_collection_stats()
    """

    # IVF_FLAT index params shared across all collections
    INDEX_PARAMS = {
        "metric_type": "COSINE",
        "index_type": "IVF_FLAT",
        "params": {"nlist": 1024},
    }

    SEARCH_PARAMS = {
        "metric_type": "COSINE",
        "params": {"nprobe": 16},
    }

    def __init__(
        self,
        host: Optional[str] = None,
        port: Optional[int] = None,
        embedding_dim: int = EMBEDDING_DIM,
    ):
        """Initialize the collection manager.

        Args:
            host: Milvus server host. Defaults to MILVUS_HOST env var or localhost.
            port: Milvus server port. Defaults to MILVUS_PORT env var or 19530.
            embedding_dim: Embedding vector dimension (384 for BGE-small-en-v1.5).
        """
        self.host = host or os.environ.get("MILVUS_HOST", "localhost")
        self.port = port or int(os.environ.get("MILVUS_PORT", "19530"))
        self.embedding_dim = embedding_dim
        self._collections: Dict[str, Collection] = {}
        self._executor = ThreadPoolExecutor(max_workers=14)

    def connect(self) -> None:
        """Connect to the Milvus server."""
        logger.info(f"Connecting to Milvus at {self.host}:{self.port}")
        connections.connect(
            alias="default",
            host=self.host,
            port=self.port,
        )
        logger.info("Connected to Milvus")

    def close(self):
        """Shut down the thread pool executor."""
        if hasattr(self, '_executor'):
            self._executor.shutdown(wait=True)

    def disconnect(self) -> None:
        """Disconnect from the Milvus server."""
        self.close()
        connections.disconnect("default")
        self._collections.clear()
        logger.info("Disconnected from Milvus")

    # ── Collection lifecycle ─────────────────────────────────────────

    def create_collection(
        self,
        name: str,
        schema: CollectionSchema,
        drop_existing: bool = False,
    ) -> Collection:
        """Create a single collection with IVF_FLAT index on the embedding field.

        Args:
            name: Collection name (must be a recognized biomarker or genomic collection).
            schema: The CollectionSchema defining the fields.
            drop_existing: If True, drop the collection first if it already exists.

        Returns:
            The pymilvus Collection object.
        """
        if drop_existing and utility.has_collection(name):
            logger.warning(f"Dropping existing collection: {name}")
            utility.drop_collection(name)

        if utility.has_collection(name):
            logger.info(f"Collection '{name}' already exists, loading reference")
            collection = Collection(name)
            self._collections[name] = collection
            return collection

        logger.info(f"Creating collection: {name}")
        collection = Collection(name=name, schema=schema)

        # Create IVF_FLAT index on the embedding field
        logger.info(f"Creating IVF_FLAT/COSINE index on '{name}.embedding'")
        collection.create_index(
            field_name="embedding",
            index_params=self.INDEX_PARAMS,
        )

        self._collections[name] = collection
        logger.info(f"Collection '{name}' created with index")
        return collection

    def create_all_collections(
        self, drop_existing: bool = False
    ) -> Dict[str, Collection]:
        """Create all 11 biomarker collections (10 domain + 1 read-only genomic).

        Args:
            drop_existing: If True, drop and recreate each collection.

        Returns:
            Dict mapping collection name to Collection object.
        """
        logger.info(f"Creating all {len(COLLECTION_SCHEMAS)} biomarker collections")
        for name, schema in COLLECTION_SCHEMAS.items():
            self.create_collection(name, schema, drop_existing=drop_existing)
        logger.info(f"All {len(COLLECTION_SCHEMAS)} collections ready")
        return dict(self._collections)

    def drop_collection(self, name: str) -> None:
        """Drop a collection by name.

        Args:
            name: The collection name to drop.
        """
        if utility.has_collection(name):
            utility.drop_collection(name)
            self._collections.pop(name, None)
            logger.info(f"Collection '{name}' dropped")
        else:
            logger.warning(f"Collection '{name}' does not exist, nothing to drop")

    def get_collection(self, name: str) -> Collection:
        """Get a collection reference, creating it if needed.

        Args:
            name: The collection name.

        Returns:
            The pymilvus Collection object.

        Raises:
            ValueError: If the name is not a recognized biomarker collection.
        """
        if name in self._collections:
            return self._collections[name]

        if utility.has_collection(name):
            collection = Collection(name)
            self._collections[name] = collection
            return collection

        if name in COLLECTION_SCHEMAS:
            return self.create_collection(name, COLLECTION_SCHEMAS[name])

        raise ValueError(
            f"Unknown collection '{name}'. "
            f"Valid collections: {list(COLLECTION_SCHEMAS.keys())}"
        )

    # ── Stats ────────────────────────────────────────────────────────

    def get_collection_stats(self) -> Dict[str, int]:
        """Get row counts for all 11 biomarker collections.

        Returns:
            Dict mapping collection name to entity count.
            Collections that do not yet exist will show 0.
        """
        stats: Dict[str, int] = {}
        for name in COLLECTION_SCHEMAS:
            if utility.has_collection(name):
                collection = Collection(name)
                stats[name] = collection.num_entities
            else:
                stats[name] = 0
        return stats

    # ── Data operations ──────────────────────────────────────────────

    def _get_output_fields(self, collection_name: str) -> List[str]:
        """Return non-embedding field names for a given collection.

        Used to build the output_fields list for search results.
        Excludes the 'embedding' field since it is large and not
        needed in result payloads.

        Args:
            collection_name: The collection to get fields for.

        Returns:
            List of field name strings (e.g. ["id", "name", "text_chunk", ...]).

        Raises:
            ValueError: If the collection_name is not recognized.
        """
        if collection_name not in COLLECTION_SCHEMAS:
            raise ValueError(
                f"Unknown collection '{collection_name}'. "
                f"Valid collections: {list(COLLECTION_SCHEMAS.keys())}"
            )

        schema = COLLECTION_SCHEMAS[collection_name]
        return [
            field.name
            for field in schema.fields
            if field.name != "embedding"
        ]

    def insert_batch(
        self,
        collection_name: str,
        records: List[Dict[str, Any]],
    ) -> int:
        """Insert a batch of records into a collection.

        Each record dict must contain all required fields for the collection
        schema, including the pre-computed 'embedding' vector.

        Args:
            collection_name: Target collection name.
            records: List of dicts with field names matching the schema.

        Returns:
            Number of records successfully inserted.
        """
        try:
            collection = self.get_collection(collection_name)
            result = collection.insert(records)
            collection.flush()
            count = result.insert_count
            logger.info(f"Inserted {count} records into {collection_name}")
            return count
        except Exception as e:
            logger.error(f"Failed to insert batch into {collection_name}: {e}")
            raise

    def search(
        self,
        collection_name: str,
        query_embedding: List[float],
        top_k: int = 5,
        filter_expr: Optional[str] = None,
        score_threshold: float = 0.0,
    ) -> List[Dict[str, Any]]:
        """Search a single collection by vector similarity.

        Args:
            collection_name: The collection to search.
            query_embedding: 384-dim query vector (BGE-small-en-v1.5).
            top_k: Maximum number of results to return.
            filter_expr: Optional Milvus boolean filter expression
                (e.g. 'gene == "MTHFR"').
            score_threshold: Minimum cosine similarity score (0.0-1.0).

        Returns:
            List of dicts with 'id', 'score', 'collection', and all output fields.
        """
        try:
            collection = self.get_collection(collection_name)
            collection.load()

            output_fields = self._get_output_fields(collection_name)

            results = collection.search(
                data=[query_embedding],
                anns_field="embedding",
                param=self.SEARCH_PARAMS,
                limit=top_k,
                output_fields=output_fields,
                expr=filter_expr,
            )

            # Convert results to list of dicts
            evidence_results: List[Dict[str, Any]] = []
            for hits in results:
                for hit in hits:
                    score = hit.score  # Cosine similarity (0-1)
                    if score < score_threshold:
                        continue

                    record: Dict[str, Any] = {
                        "id": hit.id,
                        "score": score,
                        "collection": collection_name,
                    }
                    for field_name in output_fields:
                        if field_name != "id":  # Already captured above
                            record[field_name] = hit.entity.get(field_name)

                    evidence_results.append(record)

            return evidence_results

        except Exception as e:
            logger.error(f"Search failed on {collection_name}: {e}")
            return []

    def search_all(
        self,
        query_embedding: List[float],
        top_k_per_collection: int = 5,
        filter_exprs: Optional[Dict[str, str]] = None,
        score_threshold: float = 0.0,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Search ALL biomarker collections in parallel.

        Performs vector similarity search across every collection
        concurrently using a thread pool, then merges results.

        Args:
            query_embedding: 384-dim query vector (BGE-small-en-v1.5).
            top_k_per_collection: Max results per collection.
            filter_exprs: Optional dict of collection_name -> filter expression.
                Collections not in the dict get no filter.
            score_threshold: Minimum cosine similarity score (0.0-1.0).

        Returns:
            Dict mapping collection name -> list of result dicts.
        """
        collections = list(COLLECTION_SCHEMAS.keys())
        all_results: Dict[str, List[Dict[str, Any]]] = {}

        def _search_one(name: str) -> tuple:
            expr = (filter_exprs or {}).get(name)
            return name, self.search(
                collection_name=name,
                query_embedding=query_embedding,
                top_k=top_k_per_collection,
                filter_expr=expr,
                score_threshold=score_threshold,
            )

        executor = self._executor
        futures = {
            executor.submit(_search_one, name): name
            for name in collections
        }
        for future in as_completed(futures):
            coll_name = futures[future]
            try:
                name, hits = future.result()
                all_results[name] = hits
            except Exception as e:
                logger.warning(
                    f"Search failed for collection '{coll_name}': {e}"
                )
                all_results[coll_name] = []

        total = sum(len(v) for v in all_results.values())
        logger.info(
            f"Searched {len(collections)} collections, found {total} results"
        )
        return all_results
