"""Milvus collection management for Pharmacogenomics Intelligence Agent.

Manages 15 PGx collections (14 domain-specific + 1 read-only genomic):
  - pgx_gene_reference       — Pharmacogene star allele definitions & activity scores
  - pgx_drug_guidelines      — CPIC/DPWG clinical prescribing guidelines
  - pgx_drug_interactions    — Drug-gene interaction records (PharmGKB)
  - pgx_hla_hypersensitivity — HLA-mediated adverse drug reaction screening
  - pgx_phenoconversion      — Metabolic phenoconversion via drug-drug interactions
  - pgx_dosing_algorithms    — Genotype-guided dosing algorithms & formulas
  - pgx_clinical_evidence    — Published PGx clinical evidence & outcomes
  - pgx_population_data      — Population-specific allele frequency data
  - pgx_clinical_trials      — PGx-related clinical trials
  - pgx_fda_labels           — FDA pharmacogenomic labeling information
  - pgx_drug_alternatives    — Genotype-guided therapeutic alternatives
  - pgx_patient_profiles     — Patient diplotype-phenotype profiles
  - pgx_implementation       — Clinical PGx implementation programs
  - pgx_education            — PGx educational resources & guidelines

Follows the same pymilvus pattern as:
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


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION SCHEMA DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

EMBEDDING_DIM = 384  # BGE-small-en-v1.5

# ── pgx_gene_reference ──────────────────────────────────────────────

GENE_REFERENCE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Gene-allele identifier (e.g. CYP2D6_star1)",
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
        description="Pharmacogene symbol (e.g. CYP2D6, TPMT)",
    ),
    FieldSchema(
        name="star_allele",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Star allele designation (e.g. *1, *2, *17)",
    ),
    FieldSchema(
        name="defining_variants",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Comma-separated rsIDs defining the allele",
    ),
    FieldSchema(
        name="activity_score",
        dtype=DataType.FLOAT,
        description="Enzyme activity score (e.g. 0.0, 0.5, 1.0, 2.0)",
    ),
    FieldSchema(
        name="function_status",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="no_function, decreased_function, normal_function, increased_function",
    ),
    FieldSchema(
        name="allele_frequency_global",
        dtype=DataType.FLOAT,
        description="Global allele frequency",
    ),
    FieldSchema(
        name="allele_frequency_european",
        dtype=DataType.FLOAT,
        description="European population allele frequency",
    ),
    FieldSchema(
        name="allele_frequency_african",
        dtype=DataType.FLOAT,
        description="African population allele frequency",
    ),
    FieldSchema(
        name="allele_frequency_east_asian",
        dtype=DataType.FLOAT,
        description="East Asian population allele frequency",
    ),
    FieldSchema(
        name="allele_frequency_south_asian",
        dtype=DataType.FLOAT,
        description="South Asian population allele frequency",
    ),
    FieldSchema(
        name="allele_frequency_latino",
        dtype=DataType.FLOAT,
        description="Latino population allele frequency",
    ),
    FieldSchema(
        name="pharmvar_id",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="PharmVar accession identifier",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Data source (PharmVar, CPIC, PharmGKB)",
    ),
]

GENE_REFERENCE_SCHEMA = CollectionSchema(
    fields=GENE_REFERENCE_FIELDS,
    description="Pharmacogene star allele definitions, activity scores, and population frequencies",
)

# ── pgx_drug_guidelines ─────────────────────────────────────────────

DRUG_GUIDELINES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Guideline identifier (e.g. CPIC_CYP2D6_codeine)",
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
        description="Pharmacogene symbol",
    ),
    FieldSchema(
        name="drug",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Drug name",
    ),
    FieldSchema(
        name="phenotype",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Metabolizer phenotype (e.g. poor_metabolizer, ultrarapid_metabolizer)",
    ),
    FieldSchema(
        name="guideline_body",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="CPIC, DPWG, FDA, CPNDS",
    ),
    FieldSchema(
        name="cpic_level",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="CPIC evidence level (A, A/B, B, C, D)",
    ),
    FieldSchema(
        name="recommendation",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Clinical recommendation text",
    ),
    FieldSchema(
        name="clinical_action",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Required action (e.g. avoid, dose_reduce, alternative)",
    ),
    FieldSchema(
        name="alert_level",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="high, moderate, informative",
    ),
    FieldSchema(
        name="alternative_drugs",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated alternative drug options",
    ),
    FieldSchema(
        name="dose_adjustment",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Dose modification guidance",
    ),
    FieldSchema(
        name="evidence_pmids",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated PMIDs supporting the guideline",
    ),
    FieldSchema(
        name="guideline_version",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Guideline version identifier",
    ),
    FieldSchema(
        name="last_updated",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Last update date (YYYY-MM-DD)",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

DRUG_GUIDELINES_SCHEMA = CollectionSchema(
    fields=DRUG_GUIDELINES_FIELDS,
    description="CPIC/DPWG/FDA clinical prescribing guidelines for pharmacogene-drug pairs",
)

# ── pgx_drug_interactions ───────────────────────────────────────────

DRUG_INTERACTIONS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Drug-gene interaction identifier",
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
        max_length=200,
        description="Drug name",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Pharmacogene symbol",
    ),
    FieldSchema(
        name="variant_rsid",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Specific variant rsID (e.g. rs1065852)",
    ),
    FieldSchema(
        name="interaction_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="PK, PD, efficacy, toxicity, dosage",
    ),
    FieldSchema(
        name="effect_description",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Description of the pharmacogenomic effect",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="1A, 1B, 2A, 2B, 3, 4 (PharmGKB levels)",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Clinical significance of the interaction",
    ),
    FieldSchema(
        name="pharmgkb_id",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="PharmGKB annotation identifier",
    ),
    FieldSchema(
        name="affected_phenotype",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Affected metabolizer phenotype",
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
    description="Drug-gene interaction records from PharmGKB and clinical annotations",
)

# ── pgx_hla_hypersensitivity ────────────────────────────────────────

HLA_HYPERSENSITIVITY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="HLA-drug hypersensitivity record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="hla_allele",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="HLA allele (e.g. HLA-B*57:01, HLA-A*31:01)",
    ),
    FieldSchema(
        name="drug",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Drug associated with hypersensitivity",
    ),
    FieldSchema(
        name="reaction_type",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="SJS/TEN, DRESS, HSR, hepatotoxicity",
    ),
    FieldSchema(
        name="risk_if_positive",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Risk description if HLA allele is positive",
    ),
    FieldSchema(
        name="severity",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="mild, moderate, severe, life_threatening",
    ),
    FieldSchema(
        name="cpic_level",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="CPIC evidence level (A, A/B, B)",
    ),
    FieldSchema(
        name="recommendation",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Clinical recommendation for positive carriers",
    ),
    FieldSchema(
        name="screening_mandatory",
        dtype=DataType.BOOL,
        description="Whether pre-treatment HLA screening is mandatory",
    ),
    FieldSchema(
        name="prevalence_european",
        dtype=DataType.FLOAT,
        description="HLA allele prevalence in European population",
    ),
    FieldSchema(
        name="prevalence_african",
        dtype=DataType.FLOAT,
        description="HLA allele prevalence in African population",
    ),
    FieldSchema(
        name="prevalence_east_asian",
        dtype=DataType.FLOAT,
        description="HLA allele prevalence in East Asian population",
    ),
    FieldSchema(
        name="prevalence_south_asian",
        dtype=DataType.FLOAT,
        description="HLA allele prevalence in South Asian population",
    ),
    FieldSchema(
        name="prevalence_latino",
        dtype=DataType.FLOAT,
        description="HLA allele prevalence in Latino population",
    ),
    FieldSchema(
        name="alternative_drugs",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated alternative drug options",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

HLA_HYPERSENSITIVITY_SCHEMA = CollectionSchema(
    fields=HLA_HYPERSENSITIVITY_FIELDS,
    description="HLA-mediated drug hypersensitivity screening and risk assessment",
)

# ── pgx_phenoconversion ─────────────────────────────────────────────

PHENOCONVERSION_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Phenoconversion record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="affected_enzyme",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Enzyme subject to phenoconversion (e.g. CYP2D6)",
    ),
    FieldSchema(
        name="precipitant_drug",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Drug causing the inhibition/induction",
    ),
    FieldSchema(
        name="interaction_type",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="strong_inhibitor, moderate_inhibitor, inducer",
    ),
    FieldSchema(
        name="effect_on_phenotype",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Phenotype shift (e.g. NM to PM, EM to IM)",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Clinical impact of the phenoconversion",
    ),
    FieldSchema(
        name="affected_substrate_drugs",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated substrate drugs affected",
    ),
    FieldSchema(
        name="time_to_onset",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Time to phenoconversion effect (e.g. 2-3 days)",
    ),
    FieldSchema(
        name="reversibility",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="reversible, irreversible, partially_reversible",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="strong, moderate, weak",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

PHENOCONVERSION_SCHEMA = CollectionSchema(
    fields=PHENOCONVERSION_FIELDS,
    description="Metabolic phenoconversion via drug-drug interaction-mediated enzyme inhibition/induction",
)

# ── pgx_dosing_algorithms ───────────────────────────────────────────

DOSING_ALGORITHMS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Dosing algorithm identifier",
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
        max_length=200,
        description="Drug for dosing algorithm",
    ),
    FieldSchema(
        name="genes_involved",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Comma-separated pharmacogenes used in algorithm",
    ),
    FieldSchema(
        name="algorithm_name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Algorithm name (e.g. IWPC warfarin algorithm)",
    ),
    FieldSchema(
        name="input_variables",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Required input variables (genotype, age, weight, etc.)",
    ),
    FieldSchema(
        name="formula_description",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Description of the dosing formula or decision logic",
    ),
    FieldSchema(
        name="validation_cohort",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Cohort used for algorithm validation",
    ),
    FieldSchema(
        name="accuracy_metrics",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="R-squared, MAE, or other accuracy metrics",
    ),
    FieldSchema(
        name="clinical_context",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Clinical context and applicability",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

DOSING_ALGORITHMS_SCHEMA = CollectionSchema(
    fields=DOSING_ALGORITHMS_FIELDS,
    description="Genotype-guided dosing algorithms and pharmacogenomic formulas",
)

# ── pgx_clinical_evidence ───────────────────────────────────────────

CLINICAL_EVIDENCE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Clinical evidence record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Publication or study title",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
    FieldSchema(
        name="study_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="RCT, observational, meta_analysis, case_report",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Pharmacogene studied",
    ),
    FieldSchema(
        name="drug",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Drug studied",
    ),
    FieldSchema(
        name="phenotype",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Metabolizer phenotype studied",
    ),
    FieldSchema(
        name="outcome_measure",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Primary outcome measure",
    ),
    FieldSchema(
        name="outcome_value",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Outcome value or effect size",
    ),
    FieldSchema(
        name="sample_size",
        dtype=DataType.INT64,
        description="Study sample size",
    ),
    FieldSchema(
        name="pmid",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="PubMed identifier",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Publication year",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Data source (PubMed, PharmGKB, manual)",
    ),
]

CLINICAL_EVIDENCE_SCHEMA = CollectionSchema(
    fields=CLINICAL_EVIDENCE_FIELDS,
    description="Published pharmacogenomic clinical evidence and outcomes studies",
)

# ── pgx_population_data ─────────────────────────────────────────────

POPULATION_DATA_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Population data record identifier",
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
        description="Pharmacogene symbol",
    ),
    FieldSchema(
        name="star_allele",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Star allele designation",
    ),
    FieldSchema(
        name="population",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Population or ethnicity",
    ),
    FieldSchema(
        name="allele_frequency",
        dtype=DataType.FLOAT,
        description="Allele frequency in this population",
    ),
    FieldSchema(
        name="sample_size",
        dtype=DataType.INT64,
        description="Sample size of the frequency study",
    ),
    FieldSchema(
        name="source_study",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Source study or database",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

POPULATION_DATA_SCHEMA = CollectionSchema(
    fields=POPULATION_DATA_FIELDS,
    description="Population-specific pharmacogene allele frequency data",
)

# ── pgx_clinical_trials ─────────────────────────────────────────────

CLINICAL_TRIALS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=20,
        description="NCT number (e.g. NCT03958656)",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Official trial title",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Brief summary for embedding",
    ),
    FieldSchema(
        name="nct_id",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="phase",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Trial phase",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Recruitment status",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Pharmacogene studied",
    ),
    FieldSchema(
        name="drug",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Drug studied",
    ),
    FieldSchema(
        name="enrollment",
        dtype=DataType.INT64,
        description="Target enrollment count",
    ),
    FieldSchema(
        name="start_year",
        dtype=DataType.INT64,
        description="Study start year",
    ),
    FieldSchema(
        name="outcome_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Outcome summary if available",
    ),
]

CLINICAL_TRIALS_SCHEMA = CollectionSchema(
    fields=CLINICAL_TRIALS_FIELDS,
    description="Pharmacogenomic clinical trials from ClinicalTrials.gov",
)

# ── pgx_fda_labels ──────────────────────────────────────────────────

FDA_LABELS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="FDA label record identifier",
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
        max_length=200,
        description="Drug name",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Pharmacogene referenced in label",
    ),
    FieldSchema(
        name="labeling_section",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Label section (boxed_warning, contraindication, dosage, clinical_pharmacology)",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk from FDA label",
    ),
    FieldSchema(
        name="label_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="actionable, informative, testing_required",
    ),
    FieldSchema(
        name="requirement_level",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="required, recommended, informational",
    ),
    FieldSchema(
        name="last_updated",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Last label update date (YYYY-MM-DD)",
    ),
]

FDA_LABELS_SCHEMA = CollectionSchema(
    fields=FDA_LABELS_FIELDS,
    description="FDA pharmacogenomic drug labeling information and requirements",
)

# ── pgx_drug_alternatives ───────────────────────────────────────────

DRUG_ALTERNATIVES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Drug alternative record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="primary_drug",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Primary drug requiring alternative",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Pharmacogene driving the alternative",
    ),
    FieldSchema(
        name="phenotype",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Metabolizer phenotype requiring alternative",
    ),
    FieldSchema(
        name="alternative_drug",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Recommended alternative drug",
    ),
    FieldSchema(
        name="rationale",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Clinical rationale for the alternative",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Evidence level supporting the alternative",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

DRUG_ALTERNATIVES_SCHEMA = CollectionSchema(
    fields=DRUG_ALTERNATIVES_FIELDS,
    description="Genotype-guided therapeutic alternatives for pharmacogene-drug pairs",
)

# ── pgx_patient_profiles ────────────────────────────────────────────

PATIENT_PROFILES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Patient profile record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="patient_id",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="De-identified patient identifier",
    ),
    FieldSchema(
        name="gene",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Pharmacogene symbol",
    ),
    FieldSchema(
        name="diplotype",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Diplotype (e.g. *1/*2, *1/*17)",
    ),
    FieldSchema(
        name="phenotype",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Inferred metabolizer phenotype",
    ),
    FieldSchema(
        name="activity_score",
        dtype=DataType.FLOAT,
        description="Combined activity score for the diplotype",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
]

PATIENT_PROFILES_SCHEMA = CollectionSchema(
    fields=PATIENT_PROFILES_FIELDS,
    description="Patient pharmacogenomic diplotype-to-phenotype profiles",
)

# ── pgx_implementation ──────────────────────────────────────────────

IMPLEMENTATION_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Implementation record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="institution",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Implementing institution or health system",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Implementation program title",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
    FieldSchema(
        name="workflow_type",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="preemptive, reactive, hybrid, point_of_care",
    ),
    FieldSchema(
        name="outcomes",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Reported implementation outcomes",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Implementation year",
    ),
]

IMPLEMENTATION_SCHEMA = CollectionSchema(
    fields=IMPLEMENTATION_FIELDS,
    description="Clinical pharmacogenomics implementation programs and outcomes",
)

# ── pgx_education ───────────────────────────────────────────────────

EDUCATION_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Education resource identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Educational resource title",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
    FieldSchema(
        name="audience",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Target audience (clinician, pharmacist, patient, researcher)",
    ),
    FieldSchema(
        name="topic",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Educational topic area",
    ),
    FieldSchema(
        name="complexity_level",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="introductory, intermediate, advanced",
    ),
]

EDUCATION_SCHEMA = CollectionSchema(
    fields=EDUCATION_FIELDS,
    description="Pharmacogenomics educational resources and training materials",
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
    "pgx_gene_reference": GENE_REFERENCE_SCHEMA,
    "pgx_drug_guidelines": DRUG_GUIDELINES_SCHEMA,
    "pgx_drug_interactions": DRUG_INTERACTIONS_SCHEMA,
    "pgx_hla_hypersensitivity": HLA_HYPERSENSITIVITY_SCHEMA,
    "pgx_phenoconversion": PHENOCONVERSION_SCHEMA,
    "pgx_dosing_algorithms": DOSING_ALGORITHMS_SCHEMA,
    "pgx_clinical_evidence": CLINICAL_EVIDENCE_SCHEMA,
    "pgx_population_data": POPULATION_DATA_SCHEMA,
    "pgx_clinical_trials": CLINICAL_TRIALS_SCHEMA,
    "pgx_fda_labels": FDA_LABELS_SCHEMA,
    "pgx_drug_alternatives": DRUG_ALTERNATIVES_SCHEMA,
    "pgx_patient_profiles": PATIENT_PROFILES_SCHEMA,
    "pgx_implementation": IMPLEMENTATION_SCHEMA,
    "pgx_education": EDUCATION_SCHEMA,
    "genomic_evidence": GENOMIC_EVIDENCE_SCHEMA,
}


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MANAGER
# ═══════════════════════════════════════════════════════════════════════


class PGxCollectionManager:
    """Manages 15 PGx Milvus collections (14 owned + 1 read-only genomic).

    Provides create/drop/insert/search operations across the full set of
    pharmacogenomics domain collections, following the same pymilvus patterns
    as rag-chat-pipeline/src/milvus_client.py.

    Usage:
        manager = PGxCollectionManager()
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

    def connect(self) -> None:
        """Connect to the Milvus server."""
        logger.info(f"Connecting to Milvus at {self.host}:{self.port}")
        connections.connect(
            alias="default",
            host=self.host,
            port=self.port,
        )
        logger.info("Connected to Milvus")

    def disconnect(self) -> None:
        """Disconnect from the Milvus server."""
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
            name: Collection name (must be a recognized PGx or genomic collection).
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

    def create_all_collections(self, drop_existing: bool = False) -> Dict[str, Collection]:
        """Create all 15 PGx collections (14 domain + 1 read-only genomic).

        Args:
            drop_existing: If True, drop and recreate each collection.

        Returns:
            Dict mapping collection name to Collection object.
        """
        logger.info(f"Creating all {len(COLLECTION_SCHEMAS)} PGx collections")
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
            ValueError: If the name is not a recognized PGx collection.
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
        """Get row counts for all 15 PGx collections.

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
            List of field name strings (e.g. ["id", "gene", "text_chunk", ...]).

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
                (e.g. 'gene == "CYP2D6"').
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
        """Search ALL PGx collections in parallel.

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

        with ThreadPoolExecutor(max_workers=len(collections)) as executor:
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
