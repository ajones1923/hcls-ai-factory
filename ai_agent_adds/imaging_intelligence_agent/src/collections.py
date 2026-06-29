"""Milvus collection manager for Imaging Intelligence Agent.

Manages 10 imaging-specific collections + 1 read-only genomic collection.
Mirrors the CARTCollectionManager pattern from cart_intelligence_agent.

Collections:
  - imaging_literature      — Published research papers & reviews
  - imaging_trials          — ClinicalTrials.gov AI-in-imaging records
  - imaging_findings        — Imaging finding templates & patterns
  - imaging_protocols       — Acquisition protocols & parameters
  - imaging_devices         — FDA-cleared AI/ML medical devices
  - imaging_anatomy         — Anatomical structure references
  - imaging_benchmarks      — Model performance benchmarks
  - imaging_guidelines      — Clinical practice guidelines (ACR, RSNA, etc.)
  - imaging_report_templates — Structured radiology report templates
  - imaging_datasets        — Public imaging datasets (TCIA, PhysioNet, etc.)

Read-only (created by rag-chat-pipeline):
  - genomic_evidence        — Genomic variant evidence

Follows the same pymilvus pattern as:
  rag-chat-pipeline/src/milvus_client.py (MilvusClient)
  ai_agent_adds/cart_intelligence_agent/src/collections.py (CARTCollectionManager)

Author: Adam Jones
Date: February 2026
"""

import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

from loguru import logger
from pymilvus import (
    Collection,
    CollectionSchema,
    DataType,
    FieldSchema,
    MilvusException,
    connections,
    utility,
)

from src.models import (
    AnatomyRecord,
    BenchmarkRecord,
    DatasetRecord,
    GuidelineRecord,
    ImagingDevice,
    ImagingFinding,
    ImagingLiterature,
    ImagingProtocol,
    ImagingTrial,
    ReportTemplate,
)


# ═══════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════

EMBEDDING_DIM = 384  # BGE-small-en-v1.5

INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}

SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION SCHEMA DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

# ── imaging_literature ────────────────────────────────────────────────

LITERATURE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
        description="PMID or DOI",
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
        description="Paper title",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
    FieldSchema(
        name="source_type",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="pubmed, pmc, preprint, guideline, manual",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Publication year",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality (ct, mri, xray, etc.)",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="ai_task",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="AI task type (detection, segmentation, etc.)",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Disease or indication",
    ),
    FieldSchema(
        name="keywords",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Comma-separated keywords / MeSH terms",
    ),
    FieldSchema(
        name="journal",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Journal name",
    ),
]

LITERATURE_SCHEMA = CollectionSchema(
    fields=LITERATURE_FIELDS,
    description="Medical imaging published literature and reviews",
)

# ── imaging_trials ────────────────────────────────────────────────────

TRIALS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
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
        name="sponsor",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Lead sponsor",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="ai_task",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="AI task type",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Disease or indication",
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

TRIALS_SCHEMA = CollectionSchema(
    fields=TRIALS_FIELDS,
    description="AI-in-imaging clinical trials from ClinicalTrials.gov",
)

# ── imaging_findings ──────────────────────────────────────────────────

FINDINGS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
        description="Finding record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Finding description for embedding",
    ),
    FieldSchema(
        name="finding_category",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="hemorrhage, nodule, mass, fracture, etc.",
    ),
    FieldSchema(
        name="severity",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="critical, urgent, significant, routine, normal",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="clinical_significance",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Clinical significance of the finding",
    ),
    FieldSchema(
        name="differential_diagnosis",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Differential diagnosis considerations",
    ),
    FieldSchema(
        name="recommended_followup",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Recommended follow-up actions",
    ),
    FieldSchema(
        name="measurement_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., diameter, volume, HU",
    ),
    FieldSchema(
        name="measurement_value",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., 12mm, 45 HU",
    ),
    FieldSchema(
        name="classification_system",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., Lung-RADS, BI-RADS, LI-RADS",
    ),
    FieldSchema(
        name="classification_score",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., 4A, 3",
    ),
]

FINDINGS_SCHEMA = CollectionSchema(
    fields=FINDINGS_FIELDS,
    description="Imaging finding templates and patterns",
)

# ── imaging_protocols ─────────────────────────────────────────────────

PROTOCOLS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
        description="Protocol identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Protocol description for embedding",
    ),
    FieldSchema(
        name="protocol_name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Protocol name",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="contrast_agent",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="e.g., iodinated, gadolinium, none",
    ),
    FieldSchema(
        name="slice_thickness_mm",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="e.g., 1.0, 3.0",
    ),
    FieldSchema(
        name="radiation_dose",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., CTDIvol 8 mGy",
    ),
    FieldSchema(
        name="scan_duration",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., 30 seconds",
    ),
    FieldSchema(
        name="clinical_indication",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Clinical indication for the protocol",
    ),
    FieldSchema(
        name="preprocessing_steps",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="e.g., windowing, resampling, normalization",
    ),
]

PROTOCOLS_SCHEMA = CollectionSchema(
    fields=PROTOCOLS_FIELDS,
    description="Imaging acquisition protocols and parameters",
)

# ── imaging_devices ───────────────────────────────────────────────────

DEVICES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="510(k) number or device identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Device description for embedding",
    ),
    FieldSchema(
        name="device_name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Device product name",
    ),
    FieldSchema(
        name="manufacturer",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Device manufacturer",
    ),
    FieldSchema(
        name="regulatory_status",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="510k_cleared, de_novo, pma, breakthrough, pending",
    ),
    FieldSchema(
        name="clearance_date",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="YYYY-MM-DD or YYYY-MM",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="ai_task",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="AI task type",
    ),
    FieldSchema(
        name="intended_use",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="FDA-cleared intended use statement",
    ),
    FieldSchema(
        name="performance_summary",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Performance summary from regulatory submission",
    ),
    FieldSchema(
        name="model_architecture",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Neural network architecture (e.g., 3d_unet, swin)",
    ),
]

DEVICES_SCHEMA = CollectionSchema(
    fields=DEVICES_FIELDS,
    description="FDA-cleared AI/ML medical imaging devices",
)

# ── imaging_anatomy ───────────────────────────────────────────────────

ANATOMY_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
        description="Anatomy record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Anatomy description for embedding",
    ),
    FieldSchema(
        name="structure_name",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Anatomical structure name",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="system",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., cardiovascular, respiratory, musculoskeletal",
    ),
    FieldSchema(
        name="snomed_code",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="SNOMED-CT code",
    ),
    FieldSchema(
        name="fma_id",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Foundational Model of Anatomy ID",
    ),
    FieldSchema(
        name="imaging_characteristics",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Typical imaging appearance characteristics",
    ),
    FieldSchema(
        name="common_pathologies",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Common pathologies affecting this structure",
    ),
    FieldSchema(
        name="segmentation_label_id",
        dtype=DataType.INT64,
        description="VISTA-3D segmentation label ID for Phase 2",
    ),
]

ANATOMY_SCHEMA = CollectionSchema(
    fields=ANATOMY_FIELDS,
    description="Anatomical structure references for imaging",
)

# ── imaging_benchmarks ────────────────────────────────────────────────

BENCHMARKS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
        description="Benchmark record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Benchmark description for embedding",
    ),
    FieldSchema(
        name="model_name",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="AI model name",
    ),
    FieldSchema(
        name="model_architecture",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Neural network architecture",
    ),
    FieldSchema(
        name="ai_task",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="AI task type",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="dataset_name",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Benchmark dataset name",
    ),
    FieldSchema(
        name="metric_name",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., Dice, AUC, sensitivity, specificity",
    ),
    FieldSchema(
        name="metric_value",
        dtype=DataType.FLOAT,
        description="Numeric metric value",
    ),
    FieldSchema(
        name="training_data_size",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., 10000 images",
    ),
    FieldSchema(
        name="inference_time_ms",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Inference time in milliseconds",
    ),
    FieldSchema(
        name="hardware",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="e.g., A100, DGX Spark, H100",
    ),
]

BENCHMARKS_SCHEMA = CollectionSchema(
    fields=BENCHMARKS_FIELDS,
    description="AI model performance benchmarks for medical imaging",
)

# ── imaging_guidelines ────────────────────────────────────────────────

GUIDELINES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
        description="Guideline record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Guideline description for embedding",
    ),
    FieldSchema(
        name="guideline_name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Guideline name",
    ),
    FieldSchema(
        name="organization",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., ACR, RSNA, NCCN, ESR",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Publication year",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="clinical_indication",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Clinical indication",
    ),
    FieldSchema(
        name="classification_system",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., Lung-RADS, BI-RADS, LI-RADS, TI-RADS",
    ),
    FieldSchema(
        name="key_recommendation",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Key recommendation text",
    ),
    FieldSchema(
        name="evidence_level",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="validated, emerging, exploratory",
    ),
]

GUIDELINES_SCHEMA = CollectionSchema(
    fields=GUIDELINES_FIELDS,
    description="Clinical practice guidelines for medical imaging AI",
)

# ── imaging_report_templates ──────────────────────────────────────────

REPORT_TEMPLATES_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=50,
        description="Report template identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Template description for embedding",
    ),
    FieldSchema(
        name="template_name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Report template name",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="finding_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., pulmonary_nodule, fracture, mass",
    ),
    FieldSchema(
        name="structured_fields",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Comma-separated structured report field names",
    ),
    FieldSchema(
        name="example_report",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Example structured report text",
    ),
    FieldSchema(
        name="coding_system",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., RadLex, SNOMED-CT, ICD-10",
    ),
]

REPORT_TEMPLATES_SCHEMA = CollectionSchema(
    fields=REPORT_TEMPLATES_FIELDS,
    description="Structured radiology report templates",
)

# ── imaging_datasets ──────────────────────────────────────────────────

DATASETS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Dataset identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Dataset description for embedding",
    ),
    FieldSchema(
        name="dataset_name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Dataset name",
    ),
    FieldSchema(
        name="source",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., TCIA, PhysioNet, Kaggle, MIMIC",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Imaging modality",
    ),
    FieldSchema(
        name="body_region",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Anatomical body region",
    ),
    FieldSchema(
        name="num_studies",
        dtype=DataType.INT64,
        description="Number of studies in the dataset",
    ),
    FieldSchema(
        name="num_images",
        dtype=DataType.INT64,
        description="Number of images in the dataset",
    ),
    FieldSchema(
        name="disease_labels",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Comma-separated disease labels",
    ),
    FieldSchema(
        name="annotation_type",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="e.g., bounding_box, segmentation_mask, report",
    ),
    FieldSchema(
        name="license_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., CC-BY-4.0, TCIA restricted",
    ),
    FieldSchema(
        name="download_url",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="URL to download the dataset",
    ),
]

DATASETS_SCHEMA = CollectionSchema(
    fields=DATASETS_FIELDS,
    description="Public medical imaging datasets",
)

# ── Genomic Evidence (read-only, created by rag-chat-pipeline) ───────

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

# Collections owned by the Imaging Intelligence Agent (read-write)
OWNED_COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    "imaging_literature": LITERATURE_SCHEMA,
    "imaging_trials": TRIALS_SCHEMA,
    "imaging_findings": FINDINGS_SCHEMA,
    "imaging_protocols": PROTOCOLS_SCHEMA,
    "imaging_devices": DEVICES_SCHEMA,
    "imaging_anatomy": ANATOMY_SCHEMA,
    "imaging_benchmarks": BENCHMARKS_SCHEMA,
    "imaging_guidelines": GUIDELINES_SCHEMA,
    "imaging_report_templates": REPORT_TEMPLATES_SCHEMA,
    "imaging_datasets": DATASETS_SCHEMA,
}

# All collections including read-only genomic_evidence
COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    **OWNED_COLLECTION_SCHEMAS,
    "genomic_evidence": GENOMIC_EVIDENCE_SCHEMA,
}

# Maps collection names to their Pydantic model class for validation
# genomic_evidence is None because it's read-only (no inserts from this agent)
COLLECTION_MODELS: Dict[str, type] = {
    "imaging_literature": ImagingLiterature,
    "imaging_trials": ImagingTrial,
    "imaging_findings": ImagingFinding,
    "imaging_protocols": ImagingProtocol,
    "imaging_devices": ImagingDevice,
    "imaging_anatomy": AnatomyRecord,
    "imaging_benchmarks": BenchmarkRecord,
    "imaging_guidelines": GuidelineRecord,
    "imaging_report_templates": ReportTemplate,
    "imaging_datasets": DatasetRecord,
    "genomic_evidence": None,
}

# Output fields per collection (all fields except embedding)
COLLECTION_OUTPUT_FIELDS: Dict[str, List[str]] = {
    name: [f.name for f in schema.fields if f.name != "embedding"]
    for name, schema in COLLECTION_SCHEMAS.items()
}


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MANAGER
# ═══════════════════════════════════════════════════════════════════════


class ImagingCollectionManager:
    """Manages 11 Imaging Milvus collections (10 owned + 1 read-only genomic).

    Provides create/drop/insert/search operations across the full set of
    imaging domain collections, following the same pymilvus patterns as
    rag-chat-pipeline/src/milvus_client.py and
    ai_agent_adds/cart_intelligence_agent/src/collections.py.

    Usage:
        manager = ImagingCollectionManager()
        manager.connect()
        manager.create_all_collections()
        stats = manager.get_collection_stats()
    """

    # IVF_FLAT index params shared across all collections
    INDEX_PARAMS = INDEX_PARAMS

    SEARCH_PARAMS = SEARCH_PARAMS

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
        try:
            connections.connect(
                alias="imaging",
                host=self.host,
                port=self.port,
            )
            logger.info("Connected to Milvus (alias=imaging)")
        except MilvusException as e:
            logger.error(f"Failed to connect to Milvus: {e}")
            raise

    def disconnect(self) -> None:
        """Disconnect from the Milvus server."""
        connections.disconnect("imaging")
        self._collections.clear()
        logger.info("Disconnected from Milvus (alias=imaging)")

    # ── Collection lifecycle ─────────────────────────────────────────

    def create_collection(
        self,
        name: str,
        schema: CollectionSchema,
        drop_existing: bool = False,
    ) -> Collection:
        """Create a single collection with IVF_FLAT index on the embedding field.

        Args:
            name: Collection name (must be a recognized imaging collection).
            schema: The CollectionSchema defining the fields.
            drop_existing: If True, drop the collection first if it already exists.

        Returns:
            The pymilvus Collection object.
        """
        try:
            if drop_existing and utility.has_collection(name, using="imaging"):
                logger.warning(f"Dropping existing collection: {name}")
                utility.drop_collection(name, using="imaging")

            if utility.has_collection(name, using="imaging"):
                logger.info(f"Collection '{name}' already exists, loading reference")
                collection = Collection(name, using="imaging")
                self._collections[name] = collection
                return collection

            logger.info(f"Creating collection: {name}")
            collection = Collection(name=name, schema=schema, using="imaging")

            # Create IVF_FLAT index on the embedding field
            logger.info(f"Creating IVF_FLAT/COSINE index on '{name}.embedding'")
            collection.create_index(
                field_name="embedding",
                index_params=self.INDEX_PARAMS,
            )

            self._collections[name] = collection
            logger.info(f"Collection '{name}' created with index")
            return collection

        except MilvusException as e:
            logger.error(f"Failed to create collection '{name}': {e}")
            raise

    def create_all_collections(
        self, drop_existing: bool = False
    ) -> Dict[str, Collection]:
        """Create all 10 owned imaging collections.

        Does NOT create the genomic_evidence collection (it is read-only
        and managed by the rag-chat-pipeline).

        Args:
            drop_existing: If True, drop and recreate each collection.

        Returns:
            Dict mapping collection name to Collection object.
        """
        logger.info("Creating all 10 imaging collections")
        for name, schema in OWNED_COLLECTION_SCHEMAS.items():
            self.create_collection(name, schema, drop_existing=drop_existing)
        logger.info(f"All {len(OWNED_COLLECTION_SCHEMAS)} imaging collections ready")
        return dict(self._collections)

    def drop_collection(self, name: str) -> None:
        """Drop a collection by name.

        Refuses to drop genomic_evidence (read-only, owned by rag-chat-pipeline).

        Args:
            name: The collection name to drop.
        """
        if name == "genomic_evidence":
            logger.warning(
                "Refusing to drop 'genomic_evidence' — "
                "read-only collection owned by rag-chat-pipeline"
            )
            return

        try:
            if utility.has_collection(name, using="imaging"):
                utility.drop_collection(name, using="imaging")
                self._collections.pop(name, None)
                logger.info(f"Collection '{name}' dropped")
            else:
                logger.warning(
                    f"Collection '{name}' does not exist, nothing to drop"
                )
        except MilvusException as e:
            logger.error(f"Failed to drop collection '{name}': {e}")
            raise

    def get_collection(self, name: str) -> Collection:
        """Get a collection reference, creating it if needed.

        Args:
            name: The collection name.

        Returns:
            The pymilvus Collection object.

        Raises:
            ValueError: If the name is not a recognized collection.
        """
        if name in self._collections:
            return self._collections[name]

        if utility.has_collection(name, using="imaging"):
            collection = Collection(name, using="imaging")
            self._collections[name] = collection
            return collection

        if name in OWNED_COLLECTION_SCHEMAS:
            return self.create_collection(name, OWNED_COLLECTION_SCHEMAS[name])

        raise ValueError(
            f"Unknown collection '{name}'. "
            f"Valid collections: {list(COLLECTION_SCHEMAS.keys())}"
        )

    # ── Stats ────────────────────────────────────────────────────────

    def get_collection_stats(self) -> Dict[str, int]:
        """Get row counts for all collections.

        Returns:
            Dict mapping collection name to entity count.
            Collections that do not yet exist will show 0.
        """
        stats: Dict[str, int] = {}
        for name in COLLECTION_SCHEMAS:
            try:
                if utility.has_collection(name, using="imaging"):
                    collection = Collection(name, using="imaging")
                    stats[name] = collection.num_entities
                else:
                    stats[name] = 0
            except MilvusException as e:
                logger.warning(f"Failed to get stats for '{name}': {e}")
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
            List of field name strings (e.g. ["id", "title", "text_chunk", ...]).

        Raises:
            ValueError: If the collection_name is not recognized.
        """
        if collection_name not in COLLECTION_SCHEMAS:
            raise ValueError(
                f"Unknown collection '{collection_name}'. "
                f"Valid collections: {list(COLLECTION_SCHEMAS.keys())}"
            )

        return COLLECTION_OUTPUT_FIELDS[collection_name]

    def insert_batch(
        self,
        collection_name: str,
        records: List[Dict[str, Any]],
        batch_size: int = 100,
    ) -> int:
        """Insert records into a collection in batches.

        Each record dict must contain all required fields for the collection
        schema, including the pre-computed 'embedding' vector.

        Refuses to insert into genomic_evidence (read-only).

        Args:
            collection_name: Target collection name.
            records: List of dicts with field names matching the schema.
            batch_size: Number of records per insert batch.

        Returns:
            Number of records successfully inserted.
        """
        if collection_name == "genomic_evidence":
            logger.warning(
                "Refusing to insert into 'genomic_evidence' — "
                "read-only collection owned by rag-chat-pipeline"
            )
            return 0

        if not records:
            logger.warning(f"No records to insert into {collection_name}")
            return 0

        try:
            collection = self.get_collection(collection_name)
            total_inserted = 0

            for i in range(0, len(records), batch_size):
                batch = records[i : i + batch_size]
                result = collection.insert(batch)
                total_inserted += result.insert_count
                logger.debug(
                    f"Inserted batch {i // batch_size + 1} "
                    f"({result.insert_count} records) into {collection_name}"
                )

            collection.flush()
            logger.info(
                f"Inserted {total_inserted} records into {collection_name} "
                f"({(len(records) - 1) // batch_size + 1} batches)"
            )
            return total_inserted

        except MilvusException as e:
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
                (e.g. 'modality == "ct"').
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

        except MilvusException as e:
            logger.error(f"Search failed on {collection_name}: {e}")
            return []

    def search_all(
        self,
        query_embedding: List[float],
        top_k_per_collection: int = 5,
        filter_exprs: Optional[Dict[str, str]] = None,
        score_threshold: float = 0.0,
        include_genomic: bool = False,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Search imaging collections in parallel.

        By default searches only the 10 owned imaging collections.
        Set include_genomic=True to also search the genomic_evidence
        collection.

        Args:
            query_embedding: 384-dim query vector (BGE-small-en-v1.5).
            top_k_per_collection: Max results per collection.
            filter_exprs: Optional dict of collection_name -> filter expression.
                Collections not in the dict get no filter.
            score_threshold: Minimum cosine similarity score (0.0-1.0).
            include_genomic: If True, also search genomic_evidence.

        Returns:
            Dict mapping collection name -> list of result dicts.
        """
        collections = list(OWNED_COLLECTION_SCHEMAS.keys())
        if include_genomic:
            collections.append("genomic_evidence")

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
