"""Milvus collection schemas for Cardiology Intelligence Agent.

Defines 12 domain-specific vector collections for cardiovascular medicine:
  - cardio_literature          — Published cardiology research literature
  - cardio_trials              — Cardiovascular clinical trials (ClinicalTrials.gov)
  - cardio_imaging             — Cardiac imaging protocols, findings, and criteria
  - cardio_electrophysiology   — ECG / Holter / EP / device electrophysiology
  - cardio_heart_failure       — Heart failure guidelines (HFrEF / HFmrEF / HFpEF)
  - cardio_valvular            — Valvular heart disease severity and interventions
  - cardio_prevention          — Cardiovascular prevention and risk factors
  - cardio_interventional      — Interventional / structural cardiology procedures
  - cardio_oncology            — Cardio-oncology toxicity monitoring and management
  - cardio_devices             — AI diagnostic, implantable, and wearable devices
  - cardio_guidelines          — ACC / AHA / ESC clinical practice guidelines
  - cardio_hemodynamics        — Invasive and non-invasive hemodynamic parameters

Follows the same pymilvus pattern as:
  rag-chat-pipeline/src/milvus_client.py (MilvusClient)

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

from src.models import CardioWorkflowType


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
        name: Milvus collection name (e.g. ``cardio_literature``).
        description: Human-readable description of the collection purpose.
        schema_fields: Ordered list of :class:`pymilvus.FieldSchema` objects
            defining every field in the collection (including id and embedding).
        index_params: Dict of IVF_FLAT / COSINE index parameters.
        estimated_records: Approximate number of records expected after full ingest.
        search_weight: Default relevance weight used during multi-collection search
            (0.0 – 1.0).  Higher values mean the collection contributes more to the
            final ranked result set when no workflow-specific boost is active.
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

    All 12 cardiology collections share the same embedding specification
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

# ── cardio_literature ─────────────────────────────────────────────────

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
        name="title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Publication title",
    ),
    FieldSchema(
        name="abstract",
        dtype=DataType.VARCHAR,
        max_length=8192,
        description="Full abstract text",
    ),
    FieldSchema(
        name="authors",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Comma-separated author list",
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
        name="pmid",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="PubMed identifier",
    ),
    FieldSchema(
        name="doi",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Digital Object Identifier",
    ),
    FieldSchema(
        name="mesh_terms",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Semicolon-separated MeSH terms",
    ),
    FieldSchema(
        name="study_type",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Study design (e.g. RCT, meta-analysis, cohort, case-control)",
    ),
    FieldSchema(
        name="subspecialty",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Cardiology subspecialty (e.g. heart_failure, electrophysiology, interventional)",
    ),
]

LITERATURE_CONFIG = CollectionConfig(
    name="cardio_literature",
    description="Published cardiology research literature with abstracts and MeSH terms",
    schema_fields=LITERATURE_FIELDS,
    estimated_records=3000,
    search_weight=0.10,
)

# ── cardio_trials ─────────────────────────────────────────────────────

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
        name="trial_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Official trial name or acronym (e.g. PARADIGM-HF, DAPA-HF)",
    ),
    FieldSchema(
        name="nct_id",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="ClinicalTrials.gov NCT identifier",
    ),
    FieldSchema(
        name="phase",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Trial phase (Phase 1, Phase 2, Phase 3, Phase 4)",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Target condition(s) under study",
    ),
    FieldSchema(
        name="intervention",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Drug, device, or procedural intervention",
    ),
    FieldSchema(
        name="primary_outcome",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Primary endpoint / outcome measure",
    ),
    FieldSchema(
        name="enrollment",
        dtype=DataType.INT32,
        description="Target or actual participant enrollment count",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Trial status (recruiting, completed, terminated, active)",
    ),
    FieldSchema(
        name="start_date",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Study start date (YYYY-MM-DD or YYYY-MM)",
    ),
    FieldSchema(
        name="completion_date",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Actual or estimated completion date",
    ),
    FieldSchema(
        name="sponsor",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Lead sponsor organisation",
    ),
    FieldSchema(
        name="key_findings",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Summary of key results (empty if trial still ongoing)",
    ),
]

TRIALS_CONFIG = CollectionConfig(
    name="cardio_trials",
    description="Cardiovascular clinical trials from ClinicalTrials.gov with outcomes",
    schema_fields=TRIALS_FIELDS,
    estimated_records=500,
    search_weight=0.08,
)

# ── cardio_imaging ────────────────────────────────────────────────────

IMAGING_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Imaging modality (echo, CMR, CT, nuclear, cath)",
    ),
    FieldSchema(
        name="protocol",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Imaging protocol or acquisition sequence",
    ),
    FieldSchema(
        name="finding",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Imaging finding description",
    ),
    FieldSchema(
        name="measurement_name",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Name of measured parameter (e.g. LVEF, GLS, LV mass index)",
    ),
    FieldSchema(
        name="normal_range",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Normal reference range with units",
    ),
    FieldSchema(
        name="abnormal_criteria",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Criteria indicating abnormality and severity grading",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Clinical interpretation and downstream management implications",
    ),
    FieldSchema(
        name="guideline_society",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Issuing guideline body (ASE, SCMR, SCCT, ACC/AHA)",
    ),
    FieldSchema(
        name="cross_modal_trigger",
        dtype=DataType.BOOL,
        description="Whether this finding should trigger additional imaging modality",
    ),
]

IMAGING_CONFIG = CollectionConfig(
    name="cardio_imaging",
    description="Cardiac imaging protocols, normal ranges, and abnormality criteria",
    schema_fields=IMAGING_FIELDS,
    estimated_records=200,
    search_weight=0.10,
)

# ── cardio_electrophysiology ──────────────────────────────────────────

ELECTROPHYSIOLOGY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="category",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="EP category (ECG, Holter, Device, EP)",
    ),
    FieldSchema(
        name="finding",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Electrophysiology finding or arrhythmia pattern",
    ),
    FieldSchema(
        name="criteria",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Diagnostic criteria (voltage, duration, morphology thresholds)",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Clinical significance and risk stratification",
    ),
    FieldSchema(
        name="urgency",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Clinical urgency (emergent, urgent, routine, surveillance)",
    ),
    FieldSchema(
        name="associated_conditions",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Comma-separated associated cardiac and systemic conditions",
    ),
    FieldSchema(
        name="management",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Recommended management pathway (drugs, ablation, device, observation)",
    ),
]

ELECTROPHYSIOLOGY_CONFIG = CollectionConfig(
    name="cardio_electrophysiology",
    description="ECG, Holter, EP study, and cardiac device electrophysiology data",
    schema_fields=ELECTROPHYSIOLOGY_FIELDS,
    estimated_records=150,
    search_weight=0.08,
)

# ── cardio_heart_failure ──────────────────────────────────────────────

HEART_FAILURE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="topic",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Heart failure topic heading",
    ),
    FieldSchema(
        name="hf_type",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="HF phenotype (HFrEF, HFmrEF, HFpEF)",
    ),
    FieldSchema(
        name="nyha_class",
        dtype=DataType.VARCHAR,
        max_length=8,
        description="NYHA functional class (I, II, III, IV)",
    ),
    FieldSchema(
        name="acc_stage",
        dtype=DataType.VARCHAR,
        max_length=8,
        description="ACC/AHA HF stage (A, B, C, D)",
    ),
    FieldSchema(
        name="content",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Detailed guideline content or clinical recommendation text",
    ),
    FieldSchema(
        name="drug_class",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Drug class (e.g. ARNI, SGLT2i, beta-blocker, MRA)",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Level of evidence (A, B-R, B-NR, C-LD, C-EO)",
    ),
    FieldSchema(
        name="guideline_year",
        dtype=DataType.INT32,
        description="Year of the governing guideline publication",
    ),
]

HEART_FAILURE_CONFIG = CollectionConfig(
    name="cardio_heart_failure",
    description="Heart failure management guidelines by HF type, NYHA class, and ACC stage",
    schema_fields=HEART_FAILURE_FIELDS,
    estimated_records=150,
    search_weight=0.10,
)

# ── cardio_valvular ───────────────────────────────────────────────────

VALVULAR_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="valve",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Valve name (aortic, mitral, tricuspid, pulmonic)",
    ),
    FieldSchema(
        name="pathology",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Valve pathology (stenosis, regurgitation, prolapse, endocarditis)",
    ),
    FieldSchema(
        name="severity",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Severity grade (mild, moderate, severe)",
    ),
    FieldSchema(
        name="measurement_name",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Key quantitative measurement (e.g. AVA, mean gradient, EROA)",
    ),
    FieldSchema(
        name="threshold_value",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Threshold value for severity grading with units",
    ),
    FieldSchema(
        name="intervention_type",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Intervention type (TAVR, SAVR, MitraClip, surgical repair)",
    ),
    FieldSchema(
        name="indication",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Clinical indication for intervention per guideline criteria",
    ),
    FieldSchema(
        name="guideline_reference",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Source guideline citation (e.g. 2020 ACC/AHA VHD Guideline)",
    ),
]

VALVULAR_CONFIG = CollectionConfig(
    name="cardio_valvular",
    description="Valvular heart disease severity criteria, thresholds, and intervention indications",
    schema_fields=VALVULAR_FIELDS,
    estimated_records=120,
    search_weight=0.08,
)

# ── cardio_prevention ─────────────────────────────────────────────────

PREVENTION_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="topic",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Prevention topic (e.g. lipid management, hypertension, lifestyle)",
    ),
    FieldSchema(
        name="risk_factor",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Modifiable or non-modifiable risk factor addressed",
    ),
    FieldSchema(
        name="intervention",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Preventive intervention (drug, lifestyle, screening strategy)",
    ),
    FieldSchema(
        name="target_value",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Target treatment value with units (e.g. LDL < 70 mg/dL)",
    ),
    FieldSchema(
        name="evidence_class",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="ACC/AHA class of recommendation (I, IIa, IIb, III)",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Level of evidence (A, B-R, B-NR, C-LD, C-EO)",
    ),
    FieldSchema(
        name="population",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Target patient population or risk stratum",
    ),
    FieldSchema(
        name="content",
        dtype=DataType.VARCHAR,
        max_length=4096,
        description="Full recommendation content text",
    ),
]

PREVENTION_CONFIG = CollectionConfig(
    name="cardio_prevention",
    description="Cardiovascular prevention guidelines, risk factors, and treatment targets",
    schema_fields=PREVENTION_FIELDS,
    estimated_records=150,
    search_weight=0.10,
)

# ── cardio_interventional ─────────────────────────────────────────────

INTERVENTIONAL_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="procedure_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Interventional procedure name (PCI, TAVR, PFO closure, etc.)",
    ),
    FieldSchema(
        name="indication",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Clinical indications for the procedure",
    ),
    FieldSchema(
        name="technique",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Procedural technique description and access approach",
    ),
    FieldSchema(
        name="outcomes",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Expected outcomes and success rates from landmark trials",
    ),
    FieldSchema(
        name="complications",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Known complications and their approximate incidence rates",
    ),
    FieldSchema(
        name="contraindications",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Absolute and relative contraindications",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Level of evidence (A, B-R, B-NR, C-LD, C-EO)",
    ),
]

INTERVENTIONAL_CONFIG = CollectionConfig(
    name="cardio_interventional",
    description="Interventional and structural cardiology procedures, indications, and outcomes",
    schema_fields=INTERVENTIONAL_FIELDS,
    estimated_records=100,
    search_weight=0.07,
)

# ── cardio_oncology ───────────────────────────────────────────────────

CARDIO_ONCOLOGY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="chemotherapy_agent",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Chemotherapy or immunotherapy agent name",
    ),
    FieldSchema(
        name="cardiotoxicity_type",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Type of cardiotoxicity (cardiomyopathy, arrhythmia, HTN, VTE, pericarditis)",
    ),
    FieldSchema(
        name="risk_level",
        dtype=DataType.VARCHAR,
        max_length=32,
        description="Cardiotoxicity risk level (low, moderate, high)",
    ),
    FieldSchema(
        name="monitoring_protocol",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Recommended cardiac monitoring protocol during and after treatment",
    ),
    FieldSchema(
        name="biomarker_thresholds",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Troponin / BNP / GLS thresholds triggering intervention",
    ),
    FieldSchema(
        name="management",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Management strategy for detected cardiotoxicity",
    ),
    FieldSchema(
        name="guideline_source",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Source guideline (ESC, ASCO, ACC cardio-oncology consensus)",
    ),
]

CARDIO_ONCOLOGY_CONFIG = CollectionConfig(
    name="cardio_oncology",
    description="Cardio-oncology toxicity profiles, monitoring, and management protocols",
    schema_fields=CARDIO_ONCOLOGY_FIELDS,
    estimated_records=100,
    search_weight=0.06,
)

# ── cardio_devices ────────────────────────────────────────────────────

DEVICES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="device_name",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Device or system commercial name",
    ),
    FieldSchema(
        name="device_type",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Device category (AI_diagnostic, implantable, wearable)",
    ),
    FieldSchema(
        name="manufacturer",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Device manufacturer",
    ),
    FieldSchema(
        name="fda_status",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="FDA clearance / approval status (510k, PMA, De Novo, pending)",
    ),
    FieldSchema(
        name="clinical_application",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Clinical application and use case",
    ),
    FieldSchema(
        name="evidence",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Clinical evidence summary supporting the device",
    ),
    FieldSchema(
        name="integration_notes",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="EHR / workflow integration considerations",
    ),
]

DEVICES_CONFIG = CollectionConfig(
    name="cardio_devices",
    description="Cardiac AI diagnostic, implantable, and wearable device catalogue",
    schema_fields=DEVICES_FIELDS,
    estimated_records=80,
    search_weight=0.04,
)

# ── cardio_guidelines ─────────────────────────────────────────────────

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
        name="society",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Issuing society (ACC, AHA, ESC, HRS, SCAI)",
    ),
    FieldSchema(
        name="guideline_title",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Full guideline document title",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT32,
        description="Guideline publication or update year",
    ),
    FieldSchema(
        name="recommendation",
        dtype=DataType.VARCHAR,
        max_length=2048,
        description="Individual recommendation text",
    ),
    FieldSchema(
        name="class_of_rec",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Class of recommendation (I, IIa, IIb, III)",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=16,
        description="Level of evidence (A, B-R, B-NR, C-LD, C-EO)",
    ),
    FieldSchema(
        name="condition",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Condition or disease topic addressed",
    ),
    FieldSchema(
        name="section",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Guideline section or subsection heading",
    ),
]

GUIDELINES_CONFIG = CollectionConfig(
    name="cardio_guidelines",
    description="ACC / AHA / ESC clinical practice guideline recommendations",
    schema_fields=GUIDELINES_FIELDS,
    estimated_records=150,
    search_weight=0.10,
)

# ── cardio_hemodynamics ───────────────────────────────────────────────

HEMODYNAMICS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.INT64,
        is_primary=True,
        auto_id=True,
        description="Auto-generated primary key",
    ),
    _make_embedding_field(),
    FieldSchema(
        name="parameter_name",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Hemodynamic parameter name (e.g. PCWP, CO, SVR, PA pressure)",
    ),
    FieldSchema(
        name="measurement_method",
        dtype=DataType.VARCHAR,
        max_length=256,
        description="Measurement method (Swan-Ganz, Fick, thermodilution, echo-Doppler)",
    ),
    FieldSchema(
        name="normal_range",
        dtype=DataType.VARCHAR,
        max_length=128,
        description="Normal reference range with units",
    ),
    FieldSchema(
        name="abnormal_criteria",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Criteria for mild, moderate, and severe abnormality",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=1024,
        description="Clinical significance and diagnostic implications",
    ),
    FieldSchema(
        name="calculation_formula",
        dtype=DataType.VARCHAR,
        max_length=512,
        description="Formula or derivation method (e.g. SVR = (MAP-CVP)/CO * 80)",
    ),
    FieldSchema(
        name="cathlab_context",
        dtype=DataType.VARCHAR,
        max_length=64,
        description="Context of measurement (diagnostic_cath, right_heart_cath, pre_transplant)",
    ),
]

HEMODYNAMICS_CONFIG = CollectionConfig(
    name="cardio_hemodynamics",
    description="Invasive and non-invasive hemodynamic parameters and reference ranges",
    schema_fields=HEMODYNAMICS_FIELDS,
    estimated_records=80,
    search_weight=0.06,
)


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION REGISTRY
# ═══════════════════════════════════════════════════════════════════════

ALL_COLLECTIONS: List[CollectionConfig] = [
    LITERATURE_CONFIG,
    TRIALS_CONFIG,
    IMAGING_CONFIG,
    ELECTROPHYSIOLOGY_CONFIG,
    HEART_FAILURE_CONFIG,
    VALVULAR_CONFIG,
    PREVENTION_CONFIG,
    INTERVENTIONAL_CONFIG,
    CARDIO_ONCOLOGY_CONFIG,
    DEVICES_CONFIG,
    GUIDELINES_CONFIG,
    HEMODYNAMICS_CONFIG,
]
"""Ordered list of all 12 cardiology collection configurations."""


COLLECTION_NAMES: Dict[str, str] = {
    "literature": "cardio_literature",
    "trials": "cardio_trials",
    "imaging": "cardio_imaging",
    "electrophysiology": "cardio_electrophysiology",
    "heart_failure": "cardio_heart_failure",
    "valvular": "cardio_valvular",
    "prevention": "cardio_prevention",
    "interventional": "cardio_interventional",
    "cardio_oncology": "cardio_oncology",
    "devices": "cardio_devices",
    "guidelines": "cardio_guidelines",
    "hemodynamics": "cardio_hemodynamics",
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

WORKFLOW_COLLECTION_WEIGHTS: Dict[CardioWorkflowType, Dict[str, float]] = {
    # ── Heart Failure evaluation ──────────────────────────────────────
    CardioWorkflowType.HEART_FAILURE: {
        "cardio_heart_failure": 0.25,
        "cardio_guidelines": 0.15,
        "cardio_imaging": 0.12,
        "cardio_hemodynamics": 0.10,
        "cardio_literature": 0.08,
        "cardio_trials": 0.08,
        "cardio_devices": 0.05,
        "cardio_prevention": 0.05,
        "cardio_electrophysiology": 0.04,
        "cardio_valvular": 0.04,
        "cardio_interventional": 0.02,
        "cardio_oncology": 0.02,
    },

    # ── CAD Assessment / Acute Coronary Syndrome ────────────────────
    CardioWorkflowType.CAD_ASSESSMENT: {
        "cardio_guidelines": 0.20,
        "cardio_interventional": 0.18,
        "cardio_imaging": 0.12,
        "cardio_trials": 0.10,
        "cardio_literature": 0.10,
        "cardio_electrophysiology": 0.08,
        "cardio_prevention": 0.06,
        "cardio_hemodynamics": 0.06,
        "cardio_heart_failure": 0.04,
        "cardio_devices": 0.03,
        "cardio_valvular": 0.02,
        "cardio_oncology": 0.01,
    },

    # ── Arrhythmia / electrophysiology ────────────────────────────────
    CardioWorkflowType.ARRHYTHMIA: {
        "cardio_electrophysiology": 0.25,
        "cardio_guidelines": 0.15,
        "cardio_devices": 0.12,
        "cardio_literature": 0.10,
        "cardio_trials": 0.08,
        "cardio_imaging": 0.08,
        "cardio_interventional": 0.06,
        "cardio_heart_failure": 0.05,
        "cardio_prevention": 0.04,
        "cardio_hemodynamics": 0.03,
        "cardio_valvular": 0.02,
        "cardio_oncology": 0.02,
    },

    # ── Valvular heart disease ────────────────────────────────────────
    CardioWorkflowType.VALVULAR_DISEASE: {
        "cardio_valvular": 0.25,
        "cardio_imaging": 0.18,
        "cardio_guidelines": 0.15,
        "cardio_interventional": 0.10,
        "cardio_hemodynamics": 0.08,
        "cardio_literature": 0.06,
        "cardio_trials": 0.06,
        "cardio_heart_failure": 0.04,
        "cardio_electrophysiology": 0.03,
        "cardio_devices": 0.02,
        "cardio_prevention": 0.02,
        "cardio_oncology": 0.01,
    },

    # ── Prevention / risk assessment ──────────────────────────────────
    CardioWorkflowType.PREVENTIVE_RISK: {
        "cardio_prevention": 0.25,
        "cardio_guidelines": 0.18,
        "cardio_literature": 0.12,
        "cardio_trials": 0.10,
        "cardio_imaging": 0.08,
        "cardio_devices": 0.06,
        "cardio_heart_failure": 0.05,
        "cardio_electrophysiology": 0.04,
        "cardio_interventional": 0.04,
        "cardio_hemodynamics": 0.03,
        "cardio_valvular": 0.03,
        "cardio_oncology": 0.02,
    },

    # ── Cardio-oncology consultation ──────────────────────────────────
    CardioWorkflowType.CARDIO_ONCOLOGY: {
        "cardio_oncology": 0.25,
        "cardio_imaging": 0.15,
        "cardio_guidelines": 0.12,
        "cardio_literature": 0.10,
        "cardio_heart_failure": 0.08,
        "cardio_trials": 0.08,
        "cardio_electrophysiology": 0.06,
        "cardio_devices": 0.05,
        "cardio_prevention": 0.04,
        "cardio_hemodynamics": 0.03,
        "cardio_valvular": 0.02,
        "cardio_interventional": 0.02,
    },

    # ── Stress test interpretation ──────────────────────────────────
    CardioWorkflowType.STRESS_TEST: {
        "cardio_imaging": 0.22,
        "cardio_guidelines": 0.15,
        "cardio_hemodynamics": 0.12,
        "cardio_interventional": 0.10,
        "cardio_heart_failure": 0.10,
        "cardio_literature": 0.08,
        "cardio_trials": 0.06,
        "cardio_electrophysiology": 0.05,
        "cardio_prevention": 0.04,
        "cardio_devices": 0.04,
        "cardio_valvular": 0.02,
        "cardio_oncology": 0.02,
    },

    # ── Cardiac MRI / imaging interpretation ──────────────────────────
    CardioWorkflowType.CARDIAC_MRI: {
        "cardio_imaging": 0.25,
        "cardio_guidelines": 0.12,
        "cardio_hemodynamics": 0.12,
        "cardio_heart_failure": 0.10,
        "cardio_valvular": 0.10,
        "cardio_literature": 0.08,
        "cardio_electrophysiology": 0.06,
        "cardio_trials": 0.05,
        "cardio_interventional": 0.04,
        "cardio_devices": 0.04,
        "cardio_prevention": 0.02,
        "cardio_oncology": 0.02,
    },

    # ── General cardiology query (no specific workflow) ───────────────
    CardioWorkflowType.GENERAL: {
        "cardio_literature": 0.12,
        "cardio_guidelines": 0.12,
        "cardio_imaging": 0.10,
        "cardio_heart_failure": 0.10,
        "cardio_prevention": 0.10,
        "cardio_trials": 0.08,
        "cardio_electrophysiology": 0.08,
        "cardio_valvular": 0.08,
        "cardio_interventional": 0.06,
        "cardio_hemodynamics": 0.06,
        "cardio_oncology": 0.05,
        "cardio_devices": 0.05,
    },

    # ── Acute decompensated heart failure ─────────────────────────────
    CardioWorkflowType.ACUTE_DECOMPENSATED_HF: {
        "cardio_heart_failure": 0.25,
        "cardio_hemodynamics": 0.15,
        "cardio_guidelines": 0.12,
        "cardio_imaging": 0.10,
        "cardio_devices": 0.08,
        "cardio_trials": 0.08,
        "cardio_literature": 0.06,
        "cardio_electrophysiology": 0.05,
        "cardio_prevention": 0.04,
        "cardio_valvular": 0.03,
        "cardio_interventional": 0.02,
        "cardio_oncology": 0.02,
    },

    # ── Post-MI secondary prevention ──────────────────────────────────
    CardioWorkflowType.POST_MI: {
        "cardio_guidelines": 0.18,
        "cardio_interventional": 0.18,
        "cardio_imaging": 0.12,
        "cardio_trials": 0.12,
        "cardio_prevention": 0.10,
        "cardio_heart_failure": 0.08,
        "cardio_electrophysiology": 0.06,
        "cardio_literature": 0.06,
        "cardio_hemodynamics": 0.04,
        "cardio_devices": 0.03,
        "cardio_valvular": 0.02,
        "cardio_oncology": 0.01,
    },

    # ── Myocarditis / pericarditis ────────────────────────────────────
    CardioWorkflowType.MYOCARDITIS_PERICARDITIS: {
        "cardio_imaging": 0.22,
        "cardio_guidelines": 0.15,
        "cardio_literature": 0.12,
        "cardio_heart_failure": 0.10,
        "cardio_trials": 0.10,
        "cardio_hemodynamics": 0.08,
        "cardio_electrophysiology": 0.06,
        "cardio_prevention": 0.05,
        "cardio_valvular": 0.04,
        "cardio_interventional": 0.03,
        "cardio_devices": 0.03,
        "cardio_oncology": 0.02,
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
        name: Full Milvus collection name (e.g. ``cardio_literature``)
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
    """Return a list of all 12 full Milvus collection names.

    Returns:
        Ordered list of collection name strings.
    """
    return [cfg.name for cfg in ALL_COLLECTIONS]


def get_search_weights(
    workflow: Optional[CardioWorkflowType] = None,
) -> Dict[str, float]:
    """Return collection search weights, optionally boosted for a workflow.

    When *workflow* is ``None`` (or not found in the boost table), the
    default base weights from each :class:`CollectionConfig` are returned.

    Args:
        workflow: Optional :class:`CardioWorkflowType` to apply
            workflow-specific weight boosting.

    Returns:
        Dict mapping collection name to its search weight (0.0 – 1.0).
    """
    if workflow is not None and workflow in WORKFLOW_COLLECTION_WEIGHTS:
        return dict(WORKFLOW_COLLECTION_WEIGHTS[workflow])
    return dict(_DEFAULT_SEARCH_WEIGHTS)
