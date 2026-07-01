"""Milvus collection schemas for the Neurology Intelligence Agent.

Defines 14 domain-specific vector collections for neurology clinical
decision support:
  - neuro_literature          -- Published neurology literature
  - neuro_trials              -- Clinical trials for neurological conditions
  - neuro_imaging             -- Neuroimaging findings and interpretations
  - neuro_electrophysiology   -- EEG, EMG, NCS, evoked potentials
  - neuro_degenerative        -- Neurodegenerative disease evidence
  - neuro_cerebrovascular     -- Stroke and cerebrovascular disease
  - neuro_epilepsy            -- Epilepsy syndromes and seizure data
  - neuro_oncology            -- CNS tumors and neuro-oncology
  - neuro_ms                  -- Multiple sclerosis evidence
  - neuro_movement            -- Movement disorders (PD, dystonia, etc.)
  - neuro_headache            -- Headache disorders and migraine
  - neuro_neuromuscular       -- Neuromuscular diseases
  - neuro_guidelines          -- Clinical practice guidelines
  - genomic_evidence          -- Shared genomic evidence (read-only)

Follows the same pymilvus pattern as:
  clinical_trial_intelligence_agent/src/collections.py
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

from src.models import NeuroWorkflowType


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
        name: Milvus collection name (e.g. ``neuro_literature``).
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

    All 14 neurology collections share the same embedding specification
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


def _make_id_field() -> FieldSchema:
    """Create the standard auto-increment primary key field."""
    return FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    )


# ===================================================================
# COLLECTION SCHEMA DEFINITIONS
# ===================================================================

# -- neuro_literature --------------------------------------------------

LITERATURE_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="pmid",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="PubMed ID",
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
        name="authors",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Pipe-delimited author list",
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
        name="doi",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Digital Object Identifier",
    ),
    FieldSchema(
        name="domain",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Neurology sub-domain (e.g., cerebrovascular, epilepsy)",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Evidence classification level",
    ),
    FieldSchema(
        name="study_type",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Study design (e.g., RCT, cohort, meta-analysis)",
    ),
]

LITERATURE_CONFIG = CollectionConfig(
    name="neuro_literature",
    description="Published neurology literature with abstracts, evidence levels, and domain tagging",
    schema_fields=LITERATURE_FIELDS,
    estimated_records=150000,
    search_weight=0.08,
)

# -- neuro_trials ------------------------------------------------------

TRIALS_FIELDS = [
    _make_id_field(),
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
        description="Official trial title",
    ),
    FieldSchema(
        name="summary",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Brief study summary",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Primary condition(s) under study",
    ),
    FieldSchema(
        name="intervention",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Intervention / treatment under study",
    ),
    FieldSchema(
        name="phase",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Trial phase (I, II, III, IV)",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Recruitment status",
    ),
    FieldSchema(
        name="enrollment",
        dtype=DataType.INT64,
        description="Target enrollment count",
    ),
    FieldSchema(
        name="domain",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Neurology sub-domain",
    ),
    FieldSchema(
        name="start_date",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Study start date (YYYY-MM)",
    ),
    FieldSchema(
        name="primary_outcome",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Primary outcome measure",
    ),
]

TRIALS_CONFIG = CollectionConfig(
    name="neuro_trials",
    description="Clinical trials for neurological conditions with phase, intervention, and outcome data",
    schema_fields=TRIALS_FIELDS,
    estimated_records=25000,
    search_weight=0.06,
)

# -- neuro_imaging -----------------------------------------------------

IMAGING_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Imaging modality (MRI, CT, PET, SPECT, angiography)",
    ),
    FieldSchema(
        name="sequence",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="MRI sequence or CT protocol (e.g., DWI, FLAIR, T1-contrast)",
    ),
    FieldSchema(
        name="finding",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Imaging finding description",
    ),
    FieldSchema(
        name="location",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Anatomical location of finding",
    ),
    FieldSchema(
        name="diagnosis",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Associated diagnosis or differential",
    ),
    FieldSchema(
        name="domain",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Neurology sub-domain",
    ),
    FieldSchema(
        name="urgency",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Clinical urgency (emergent, urgent, routine)",
    ),
    FieldSchema(
        name="pattern",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Characteristic imaging pattern (e.g., ring-enhancing, DWI restriction)",
    ),
    FieldSchema(
        name="reference",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Source reference or guideline",
    ),
]

IMAGING_CONFIG = CollectionConfig(
    name="neuro_imaging",
    description="Neuroimaging findings, patterns, and differential diagnoses across MRI/CT/PET modalities",
    schema_fields=IMAGING_FIELDS,
    estimated_records=50000,
    search_weight=0.09,
)

# -- neuro_electrophysiology -------------------------------------------

ELECTROPHYSIOLOGY_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="test_type",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Test type (EEG, EMG, NCS, VEP, SSEP, BAEP)",
    ),
    FieldSchema(
        name="finding",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Electrophysiology finding description",
    ),
    FieldSchema(
        name="pattern",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Characteristic pattern (e.g., 3Hz spike-wave, fibrillation potentials)",
    ),
    FieldSchema(
        name="lateralization",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Lateralization (left, right, bilateral, generalized)",
    ),
    FieldSchema(
        name="localization",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Anatomical localization (e.g., temporal, frontal, C5-T1 myotome)",
    ),
    FieldSchema(
        name="clinical_correlation",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Clinical correlation and significance",
    ),
    FieldSchema(
        name="domain",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Neurology sub-domain",
    ),
    FieldSchema(
        name="normal_values",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Normal reference values for comparison",
    ),
]

ELECTROPHYSIOLOGY_CONFIG = CollectionConfig(
    name="neuro_electrophysiology",
    description="EEG, EMG, NCS, and evoked potential findings with patterns and clinical correlations",
    schema_fields=ELECTROPHYSIOLOGY_FIELDS,
    estimated_records=30000,
    search_weight=0.07,
)

# -- neuro_degenerative ------------------------------------------------

DEGENERATIVE_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Neurodegenerative disease name",
    ),
    FieldSchema(
        name="subtype",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Disease subtype or variant",
    ),
    FieldSchema(
        name="diagnostic_criteria",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Established diagnostic criteria summary",
    ),
    FieldSchema(
        name="biomarkers",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Relevant biomarkers (CSF, blood, imaging)",
    ),
    FieldSchema(
        name="genetics",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Known genetic associations",
    ),
    FieldSchema(
        name="staging",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Disease staging system and criteria",
    ),
    FieldSchema(
        name="treatments",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Current treatment options (approved and investigational)",
    ),
    FieldSchema(
        name="prognosis",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Natural history and prognostic factors",
    ),
    FieldSchema(
        name="prevalence",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Prevalence or incidence estimate",
    ),
]

DEGENERATIVE_CONFIG = CollectionConfig(
    name="neuro_degenerative",
    description="Neurodegenerative diseases with diagnostic criteria, biomarkers, genetics, and staging",
    schema_fields=DEGENERATIVE_FIELDS,
    estimated_records=15000,
    search_weight=0.09,
)

# -- neuro_cerebrovascular ---------------------------------------------

CEREBROVASCULAR_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Cerebrovascular condition (ischemic stroke, ICH, SAH, TIA, etc.)",
    ),
    FieldSchema(
        name="subtype",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Stroke subtype or mechanism (e.g., large artery, cardioembolic)",
    ),
    FieldSchema(
        name="presentation",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Clinical presentation and symptoms",
    ),
    FieldSchema(
        name="workup",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Recommended diagnostic workup",
    ),
    FieldSchema(
        name="treatment_acute",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Acute treatment protocols (tPA, thrombectomy, BP management)",
    ),
    FieldSchema(
        name="treatment_secondary",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Secondary prevention strategies",
    ),
    FieldSchema(
        name="time_window",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Treatment time windows and eligibility criteria",
    ),
    FieldSchema(
        name="scoring_scales",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Relevant scoring scales (NIHSS, ASPECTS, mRS)",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Evidence classification level",
    ),
]

CEREBROVASCULAR_CONFIG = CollectionConfig(
    name="neuro_cerebrovascular",
    description="Stroke and cerebrovascular disease with acute management, time windows, and secondary prevention",
    schema_fields=CEREBROVASCULAR_FIELDS,
    estimated_records=20000,
    search_weight=0.09,
)

# -- neuro_epilepsy ----------------------------------------------------

EPILEPSY_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="syndrome",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Epilepsy syndrome (e.g., TLE, JME, Dravet)",
    ),
    FieldSchema(
        name="seizure_types",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Associated seizure types (ILAE classification)",
    ),
    FieldSchema(
        name="eeg_pattern",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Characteristic EEG pattern",
    ),
    FieldSchema(
        name="imaging_findings",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Typical neuroimaging findings",
    ),
    FieldSchema(
        name="genetics",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Known genetic associations (SCN1A, etc.)",
    ),
    FieldSchema(
        name="first_line_aed",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="First-line anti-epileptic drug recommendations",
    ),
    FieldSchema(
        name="surgical_options",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Surgical treatment options if applicable",
    ),
    FieldSchema(
        name="prognosis",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Long-term prognosis and seizure freedom rates",
    ),
    FieldSchema(
        name="age_onset",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Typical age of onset range",
    ),
]

EPILEPSY_CONFIG = CollectionConfig(
    name="neuro_epilepsy",
    description="Epilepsy syndromes with seizure types, EEG patterns, genetic associations, and treatment",
    schema_fields=EPILEPSY_FIELDS,
    estimated_records=12000,
    search_weight=0.08,
)

# -- neuro_oncology ----------------------------------------------------

ONCOLOGY_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="tumor_type",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="WHO 2021 CNS tumor classification",
    ),
    FieldSchema(
        name="who_grade",
        dtype=DataType.VARCHAR,
        max_length=8,
        description="WHO tumor grade (1-4)",
    ),
    FieldSchema(
        name="molecular_profile",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Key molecular markers (IDH, MGMT, 1p19q, H3K27M, etc.)",
    ),
    FieldSchema(
        name="location",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Typical tumor location(s)",
    ),
    FieldSchema(
        name="imaging_features",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Characteristic imaging features",
    ),
    FieldSchema(
        name="treatment_protocol",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Standard treatment protocol (surgery, RT, chemotherapy)",
    ),
    FieldSchema(
        name="prognosis",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Survival statistics and prognostic factors",
    ),
    FieldSchema(
        name="epidemiology",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Incidence and demographics",
    ),
]

ONCOLOGY_CONFIG = CollectionConfig(
    name="neuro_oncology",
    description="CNS tumors with WHO 2021 classification, molecular profiles, and treatment protocols",
    schema_fields=ONCOLOGY_FIELDS,
    estimated_records=8000,
    search_weight=0.06,
)

# -- neuro_ms ----------------------------------------------------------

MS_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="phenotype",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="MS phenotype (CIS, RRMS, SPMS, PPMS)",
    ),
    FieldSchema(
        name="topic",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Clinical topic (diagnosis, DMT selection, relapse management, etc.)",
    ),
    FieldSchema(
        name="content",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Evidence summary or guideline content",
    ),
    FieldSchema(
        name="dmt_name",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Disease-modifying therapy name if applicable",
    ),
    FieldSchema(
        name="dmt_category",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="DMT efficacy category (platform, moderate, high)",
    ),
    FieldSchema(
        name="monitoring",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Monitoring requirements and safety considerations",
    ),
    FieldSchema(
        name="mri_criteria",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="MRI-based diagnostic or monitoring criteria",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Evidence classification level",
    ),
]

MS_CONFIG = CollectionConfig(
    name="neuro_ms",
    description="Multiple sclerosis evidence with phenotype-specific DMT guidance and monitoring criteria",
    schema_fields=MS_FIELDS,
    estimated_records=10000,
    search_weight=0.07,
)

# -- neuro_movement ----------------------------------------------------

MOVEMENT_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="disorder",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Movement disorder name (PD, ET, dystonia, Huntington's, etc.)",
    ),
    FieldSchema(
        name="subtype",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Disease subtype or variant",
    ),
    FieldSchema(
        name="motor_features",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Motor features and phenomenology",
    ),
    FieldSchema(
        name="non_motor_features",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Non-motor symptoms and features",
    ),
    FieldSchema(
        name="diagnostic_criteria",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Diagnostic criteria summary",
    ),
    FieldSchema(
        name="genetics",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Genetic associations (LRRK2, GBA, SNCA, HTT, etc.)",
    ),
    FieldSchema(
        name="treatments",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Pharmacological and surgical treatments",
    ),
    FieldSchema(
        name="scales",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Assessment scales (MDS-UPDRS, Hoehn & Yahr, etc.)",
    ),
]

MOVEMENT_CONFIG = CollectionConfig(
    name="neuro_movement",
    description="Movement disorders with motor/non-motor features, genetics, and treatment options",
    schema_fields=MOVEMENT_FIELDS,
    estimated_records=12000,
    search_weight=0.07,
)

# -- neuro_headache ----------------------------------------------------

HEADACHE_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="headache_type",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="ICHD-3 headache classification",
    ),
    FieldSchema(
        name="diagnostic_criteria",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="ICHD-3 diagnostic criteria",
    ),
    FieldSchema(
        name="features",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Clinical features and characteristics",
    ),
    FieldSchema(
        name="acute_treatment",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Acute/abortive treatment options",
    ),
    FieldSchema(
        name="preventive_treatment",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Preventive treatment options",
    ),
    FieldSchema(
        name="red_flags",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Red flag symptoms requiring urgent workup",
    ),
    FieldSchema(
        name="differential",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Key differential diagnoses",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Evidence classification level",
    ),
]

HEADACHE_CONFIG = CollectionConfig(
    name="neuro_headache",
    description="Headache disorders with ICHD-3 criteria, treatment protocols, and red flags",
    schema_fields=HEADACHE_FIELDS,
    estimated_records=8000,
    search_weight=0.06,
)

# -- neuro_neuromuscular -----------------------------------------------

NEUROMUSCULAR_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Neuromuscular disease (ALS, MG, GBS, CIDP, myopathy, etc.)",
    ),
    FieldSchema(
        name="category",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Category (motor neuron, NMJ, neuropathy, myopathy)",
    ),
    FieldSchema(
        name="clinical_features",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Clinical features and presentation",
    ),
    FieldSchema(
        name="emg_pattern",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Characteristic EMG/NCS pattern",
    ),
    FieldSchema(
        name="antibodies",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Relevant antibodies (AChR, MuSK, LRP4, ganglioside, etc.)",
    ),
    FieldSchema(
        name="genetics",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Genetic associations (SOD1, SMN1, DMD, etc.)",
    ),
    FieldSchema(
        name="treatments",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Treatment options including disease-modifying therapies",
    ),
    FieldSchema(
        name="prognosis",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Prognosis and natural history",
    ),
    FieldSchema(
        name="diagnostic_criteria",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Diagnostic criteria (e.g., El Escorial for ALS)",
    ),
]

NEUROMUSCULAR_CONFIG = CollectionConfig(
    name="neuro_neuromuscular",
    description="Neuromuscular diseases with EMG patterns, antibodies, genetics, and diagnostic criteria",
    schema_fields=NEUROMUSCULAR_FIELDS,
    estimated_records=10000,
    search_weight=0.06,
)

# -- neuro_guidelines --------------------------------------------------

GUIDELINES_FIELDS = [
    _make_id_field(),
    _make_embedding_field(),
    FieldSchema(
        name="guideline_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Guideline identifier or code",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Guideline title",
    ),
    FieldSchema(
        name="organization",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Issuing organization (AAN, AHA/ASA, ILAE, etc.)",
    ),
    FieldSchema(
        name="recommendation",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Recommendation text",
    ),
    FieldSchema(
        name="guideline_class",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Recommendation class (I, IIa, IIb, III)",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Evidence level (A, B-R, B-NR, C-LD, C-EO)",
    ),
    FieldSchema(
        name="domain",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Neurology sub-domain",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Publication year",
    ),
    FieldSchema(
        name="doi",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Digital Object Identifier",
    ),
]

GUIDELINES_CONFIG = CollectionConfig(
    name="neuro_guidelines",
    description="Clinical practice guidelines from AAN/AHA/ILAE with recommendation class and evidence level",
    schema_fields=GUIDELINES_FIELDS,
    estimated_records=5000,
    search_weight=0.07,
)

# -- genomic_evidence (shared) -----------------------------------------

GENOMIC_FIELDS = [
    _make_id_field(),
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
        description="Variant description (HGVS notation)",
    ),
    FieldSchema(
        name="classification",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ACMG classification (pathogenic, likely pathogenic, VUS, etc.)",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Associated condition or phenotype",
    ),
    FieldSchema(
        name="evidence_summary",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Evidence summary for variant classification",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Data source (ClinVar, gnomAD, etc.)",
    ),
    FieldSchema(
        name="allele_frequency",
        dtype=DataType.FLOAT,
        description="Population allele frequency",
    ),
]

GENOMIC_CONFIG = CollectionConfig(
    name="genomic_evidence",
    description="Shared genomic evidence collection with variant classifications and gene-disease associations",
    schema_fields=GENOMIC_FIELDS,
    estimated_records=500000,
    search_weight=0.05,
)


# ===================================================================
# ALL COLLECTIONS LIST
# ===================================================================

ALL_COLLECTIONS: List[CollectionConfig] = [
    LITERATURE_CONFIG,
    TRIALS_CONFIG,
    IMAGING_CONFIG,
    ELECTROPHYSIOLOGY_CONFIG,
    DEGENERATIVE_CONFIG,
    CEREBROVASCULAR_CONFIG,
    EPILEPSY_CONFIG,
    ONCOLOGY_CONFIG,
    MS_CONFIG,
    MOVEMENT_CONFIG,
    HEADACHE_CONFIG,
    NEUROMUSCULAR_CONFIG,
    GUIDELINES_CONFIG,
    GENOMIC_CONFIG,
]
"""Ordered list of all 14 neurology collection configurations."""


COLLECTION_NAMES: Dict[str, str] = {
    "literature": "neuro_literature",
    "trials": "neuro_trials",
    "imaging": "neuro_imaging",
    "electrophysiology": "neuro_electrophysiology",
    "degenerative": "neuro_degenerative",
    "cerebrovascular": "neuro_cerebrovascular",
    "epilepsy": "neuro_epilepsy",
    "oncology": "neuro_oncology",
    "ms": "neuro_ms",
    "movement": "neuro_movement",
    "headache": "neuro_headache",
    "neuromuscular": "neuro_neuromuscular",
    "guidelines": "neuro_guidelines",
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

WORKFLOW_COLLECTION_WEIGHTS: Dict[NeuroWorkflowType, Dict[str, float]] = {

    # -- Acute Stroke --------------------------------------------------
    NeuroWorkflowType.ACUTE_STROKE: {
        "neuro_cerebrovascular": 0.25,
        "neuro_imaging": 0.18,
        "neuro_guidelines": 0.12,
        "neuro_literature": 0.10,
        "neuro_trials": 0.06,
        "neuro_electrophysiology": 0.04,
        "neuro_degenerative": 0.03,
        "neuro_epilepsy": 0.03,
        "neuro_oncology": 0.02,
        "neuro_ms": 0.02,
        "neuro_movement": 0.02,
        "neuro_headache": 0.03,
        "neuro_neuromuscular": 0.02,
        "genomic_evidence": 0.08,
    },

    # -- Dementia Evaluation -------------------------------------------
    NeuroWorkflowType.DEMENTIA_EVALUATION: {
        "neuro_degenerative": 0.25,
        "neuro_imaging": 0.15,
        "neuro_literature": 0.10,
        "neuro_guidelines": 0.10,
        "genomic_evidence": 0.08,
        "neuro_electrophysiology": 0.06,
        "neuro_trials": 0.06,
        "neuro_movement": 0.05,
        "neuro_cerebrovascular": 0.04,
        "neuro_ms": 0.03,
        "neuro_neuromuscular": 0.02,
        "neuro_epilepsy": 0.02,
        "neuro_headache": 0.02,
        "neuro_oncology": 0.02,
    },

    # -- Epilepsy Focus ------------------------------------------------
    NeuroWorkflowType.EPILEPSY_FOCUS: {
        "neuro_epilepsy": 0.25,
        "neuro_electrophysiology": 0.20,
        "neuro_imaging": 0.15,
        "neuro_guidelines": 0.08,
        "neuro_literature": 0.07,
        "genomic_evidence": 0.05,
        "neuro_trials": 0.05,
        "neuro_degenerative": 0.03,
        "neuro_oncology": 0.03,
        "neuro_cerebrovascular": 0.02,
        "neuro_movement": 0.02,
        "neuro_neuromuscular": 0.02,
        "neuro_ms": 0.02,
        "neuro_headache": 0.01,
    },

    # -- Brain Tumor ---------------------------------------------------
    NeuroWorkflowType.BRAIN_TUMOR: {
        "neuro_oncology": 0.25,
        "neuro_imaging": 0.18,
        "neuro_guidelines": 0.10,
        "neuro_literature": 0.10,
        "neuro_trials": 0.08,
        "genomic_evidence": 0.08,
        "neuro_epilepsy": 0.04,
        "neuro_cerebrovascular": 0.03,
        "neuro_degenerative": 0.03,
        "neuro_electrophysiology": 0.03,
        "neuro_headache": 0.03,
        "neuro_neuromuscular": 0.02,
        "neuro_ms": 0.02,
        "neuro_movement": 0.01,
    },

    # -- MS Monitoring -------------------------------------------------
    NeuroWorkflowType.MS_MONITORING: {
        "neuro_ms": 0.28,
        "neuro_imaging": 0.15,
        "neuro_guidelines": 0.12,
        "neuro_literature": 0.10,
        "neuro_trials": 0.08,
        "neuro_electrophysiology": 0.05,
        "genomic_evidence": 0.04,
        "neuro_degenerative": 0.04,
        "neuro_neuromuscular": 0.03,
        "neuro_cerebrovascular": 0.03,
        "neuro_epilepsy": 0.02,
        "neuro_movement": 0.02,
        "neuro_headache": 0.02,
        "neuro_oncology": 0.02,
    },

    # -- Parkinson's Assessment ----------------------------------------
    NeuroWorkflowType.PARKINSONS_ASSESSMENT: {
        "neuro_movement": 0.25,
        "neuro_degenerative": 0.18,
        "neuro_imaging": 0.12,
        "neuro_guidelines": 0.10,
        "neuro_literature": 0.08,
        "genomic_evidence": 0.06,
        "neuro_trials": 0.05,
        "neuro_electrophysiology": 0.04,
        "neuro_neuromuscular": 0.03,
        "neuro_cerebrovascular": 0.03,
        "neuro_epilepsy": 0.02,
        "neuro_ms": 0.02,
        "neuro_headache": 0.01,
        "neuro_oncology": 0.01,
    },

    # -- Headache Classification ---------------------------------------
    NeuroWorkflowType.HEADACHE_CLASSIFICATION: {
        "neuro_headache": 0.30,
        "neuro_guidelines": 0.15,
        "neuro_imaging": 0.12,
        "neuro_literature": 0.10,
        "neuro_trials": 0.06,
        "neuro_cerebrovascular": 0.06,
        "neuro_oncology": 0.04,
        "neuro_electrophysiology": 0.03,
        "neuro_ms": 0.03,
        "neuro_degenerative": 0.03,
        "neuro_epilepsy": 0.02,
        "neuro_neuromuscular": 0.02,
        "neuro_movement": 0.02,
        "genomic_evidence": 0.02,
    },

    # -- Neuromuscular Evaluation --------------------------------------
    NeuroWorkflowType.NEUROMUSCULAR_EVALUATION: {
        "neuro_neuromuscular": 0.28,
        "neuro_electrophysiology": 0.18,
        "neuro_guidelines": 0.10,
        "neuro_literature": 0.10,
        "genomic_evidence": 0.08,
        "neuro_trials": 0.06,
        "neuro_imaging": 0.05,
        "neuro_degenerative": 0.04,
        "neuro_movement": 0.03,
        "neuro_cerebrovascular": 0.02,
        "neuro_ms": 0.02,
        "neuro_epilepsy": 0.02,
        "neuro_headache": 0.01,
        "neuro_oncology": 0.01,
    },

    # -- General (no specific workflow) --------------------------------
    NeuroWorkflowType.GENERAL: {
        "neuro_literature": 0.08,
        "neuro_trials": 0.06,
        "neuro_imaging": 0.09,
        "neuro_electrophysiology": 0.07,
        "neuro_degenerative": 0.09,
        "neuro_cerebrovascular": 0.09,
        "neuro_epilepsy": 0.08,
        "neuro_oncology": 0.06,
        "neuro_ms": 0.07,
        "neuro_movement": 0.07,
        "neuro_headache": 0.06,
        "neuro_neuromuscular": 0.06,
        "neuro_guidelines": 0.07,
        "genomic_evidence": 0.05,
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
        name: Full Milvus collection name (e.g. ``neuro_literature``)
            **or** a short alias (e.g. ``literature``).

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
    workflow: Optional[NeuroWorkflowType] = None,
) -> Dict[str, float]:
    """Return collection search weights, optionally boosted for a workflow.

    When *workflow* is ``None`` (or not found in the boost table), the
    default base weights from each :class:`CollectionConfig` are returned.

    Args:
        workflow: Optional :class:`NeuroWorkflowType` to apply
            workflow-specific weight boosting.

    Returns:
        Dict mapping collection name to its search weight (0.0 - 1.0).
    """
    if workflow is not None and workflow in WORKFLOW_COLLECTION_WEIGHTS:
        return dict(WORKFLOW_COLLECTION_WEIGHTS[workflow])
    return dict(_DEFAULT_SEARCH_WEIGHTS)
