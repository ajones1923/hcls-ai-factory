# Precision Autoimmune Intelligence Agent -- Project Bible

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0
**License:** Apache 2.0

> Complete implementation reference. Import this document as context for Claude Code sessions working on the Precision Autoimmune Intelligence Agent.

---

## Table of Contents

1. [Project Overview & Goals](#1-project-overview-goals)
2. [DGX Spark Hardware Reference](#2-dgx-spark-hardware-reference)
3. [Repository Layout](#3-repository-layout)
4. [Docker Compose Services](#4-docker-compose-services)
5. [Milvus Collection Schemas](#5-milvus-collection-schemas)
6. [Pydantic Data Models](#6-pydantic-data-models)
7. [Configuration Reference](#7-configuration-reference)
8. [Embedding Strategy](#8-embedding-strategy)
9. [AutoantibodyInterpreter](#9-autoantibodyinterpreter)
10. [HLAAssociationAnalyzer](#10-hlaassociationanalyzer)
11. [DiseaseActivityScorer](#11-diseaseactivityscorer)
12. [FlarePredictor](#12-flarepredioctor)
13. [BiologicTherapyAdvisor](#13-biologictherapyadvisor)
14. [DiagnosticOdysseyAnalyzer](#14-diagnosticodysseyanalyzer)
15. [Multi-Collection RAG Engine](#15-multi-collection-rag-engine)
16. [Agent Orchestrator](#16-agent-orchestrator)
17. [Export Pipeline](#17-export-pipeline)
18. [FastAPI REST Server](#18-fastapi-rest-server)
19. [Streamlit UI](#19-streamlit-ui)
20. [Demo Patient Data](#20-demo-patient-data)
21. [Cross-Agent Integration](#21-cross-agent-integration)
22. [Monitoring & Metrics](#22-monitoring-metrics)
23. [Testing Strategy](#23-testing-strategy)
24. [HCLS AI Factory Integration](#24-hcls-ai-factory-integration)
25. [Implementation Sequence](#25-implementation-sequence)

---

## 1. Project Overview & Goals

### What This Agent Does

The Precision Autoimmune Intelligence Agent transforms fragmented autoimmune clinical data -- scattered across years of multi-specialist visits, lab reports, and imaging studies -- into unified, actionable intelligence. It integrates:

- **Autoantibody panels** (24 autoantibodies mapped to 13+ diseases with sensitivity/specificity)
- **HLA typing** (22 alleles with disease odds ratios and PubMed references)
- **Disease activity scoring** (20 validated scoring systems across 13 diseases)
- **Flare prediction** (13 disease-specific biomarker patterns)
- **Biologic therapy selection** (22 therapies with pharmacogenomic considerations)
- **Diagnostic odyssey analysis** (timeline reconstruction, overlap syndrome detection, classification criteria evaluation)
- **Clinical document ingestion** (PDF chunking, embedding, and RAG-powered search)

### Key Results

| Metric | Value |
|---|---|
| Unit tests passing | 455 |
| Milvus collections | 14 (13 owned + 1 read-only) |
| Autoimmune diseases | 13 in AutoimmuneDisease enum |
| Biologic therapies | 22 with PGx considerations |
| Autoantibodies mapped | 24 to disease associations |
| HLA alleles | 22 with disease odds ratios |
| Disease activity scores | 20 scoring systems |
| Flare patterns | 13 disease-specific |
| Classification criteria | 10 ACR/EULAR sets |
| Overlap syndromes | 9 cross-disease patterns |
| Demo patients | 9 with full clinical PDFs |
| Python files | 44 |
| Lines of code | ~20,000 |
| Knowledge version | v2.0.0 |

### Pipeline Pattern

```
Patient Data (Autoantibodies + HLA + Biomarkers + Clinical PDFs)
    |
    ├──> AutoantibodyInterpreter (24 antibodies -> disease associations)
    ├──> HLAAssociationAnalyzer (22 alleles -> odds ratios)
    ├──> DiseaseActivityScorer (20 scoring systems)
    ├──> FlarePredictor (13 disease patterns)
    ├──> BiologicTherapyAdvisor (22 therapies + PGx)
    └──> DiagnosticOdysseyAnalyzer (criteria + overlap + differential)
    |
    v
RAG Engine (14-collection search + Claude LLM)
    |
    v
Export (FHIR R4 | PDF | Markdown)
```

### VAST AI OS Integration

| VAST AI OS Component | Autoimmune Agent Role |
|---|---|
| **DataStore** | Clinical PDFs, knowledge base JSON (HLA, antibodies, biologics, flare patterns) |
| **DataEngine** | document_processor.py: PDF -> chunks -> BGE-small -> Milvus |
| **DataBase** | 14 Milvus collections + 9 demo patients |
| **InsightEngine** | 6 clinical engines + multi-collection RAG |
| **AgentEngine** | AutoimmuneAgent + Streamlit + FastAPI |

---

## 2. DGX Spark Hardware Reference

### Specifications

| Component | Specification |
|---|---|
| GPU | NVIDIA GB10 Grace Blackwell Superchip |
| GPU Memory | Unified 128 GB LPDDR5X |
| CPU | 10-core Arm Cortex |
| Storage | Up to 4 TB NVMe SSD |
| Networking | 10 GbE, WiFi 7, Bluetooth 5.4 |
| OS | Ubuntu-based Linux for DGX |
| CUDA | 12.x |
| Power | Desktop form factor, standard AC |
| Price | $3,999 |

### Why DGX Spark

- Unified memory eliminates CPU-GPU data transfers for embedding generation
- Desktop deployment enables on-premise HIPAA-compliant processing
- Sufficient for all 6 clinical engines + Milvus + Streamlit simultaneously
- BGE-small inference runs entirely on GPU with <5ms latency
- PDF ingestion with PyPDF2 + embedding pipeline runs efficiently

---

## 3. Repository Layout

```
precision_autoimmune_agent/
├── src/
│   ├── models.py                    # Pydantic models (AutoimmunePatientProfile, AutoantibodyPanel, etc.)
│   ├── collections.py               # 14 Milvus collection schemas + AutoimmuneCollectionManager
│   ├── rag_engine.py                # AutoimmuneRAGEngine (parallel search + Claude streaming)
│   ├── agent.py                     # AutoimmuneAgent orchestrator (6 analysis engines)
│   ├── knowledge.py                 # Knowledge v2.0.0 (HLA, antibodies, biologics, flare patterns)
│   ├── diagnostic_engine.py         # Classification criteria, overlap syndromes, differential diagnosis
│   ├── document_processor.py        # PDF ingestion, chunking, embedding (~435 lines)
│   ├── timeline_builder.py          # Diagnostic odyssey timeline construction (~251 lines)
│   ├── export.py                    # FHIR R4, PDF, Markdown export (~389 lines)
│   └── __init__.py
├── app/
│   └── autoimmune_ui.py             # Streamlit UI, 10 tabs (~1,126 lines)
├── api/
│   ├── main.py                      # FastAPI REST server, 14 endpoints (~583 lines)
│   └── routes/
│       └── __init__.py
├── config/
│   └── settings.py                  # AutoimmuneSettings (AUTO_ env prefix, ~243 lines)
├── data/
│   ├── reference/                   # Reference data files
│   ├── cache/                       # Embedding cache
│   └── events/                      # Cross-agent event store
├── demo_data/
│   ├── sarah_mitchell/              # SLE patient (34F) -- 27+ clinical PDFs
│   ├── maya_rodriguez/              # POTS/hEDS/MCAS patient (28F) -- clinical PDFs
│   ├── linda_chen/                  # Sjogren's patient (45F) -- clinical PDFs
│   ├── david_park/                  # AS patient (45M) -- clinical PDFs
│   ├── rachel_thompson/             # MCTD patient (38F) -- clinical PDFs
│   ├── emma_williams/               # MS (RRMS) patient (34F) -- clinical PDFs
│   ├── james_cooper/                # T1D + Celiac (19M) -- clinical PDFs
│   ├── karen_foster/                # SSc (dcSSc) patient (48F) -- clinical PDFs
│   └── michael_torres/              # Graves' Disease patient (41M) -- clinical PDFs
├── scripts/
│   ├── setup_collections.py         # Create and seed Milvus collections
│   ├── generate_demo_patients.py    # Generate demo patient PDF records
│   ├── pdf_engine.py                # PDF generation engine
│   ├── patient_sarah.py             # Sarah Mitchell patient data
│   ├── patient_emma.py              # Emma Williams patient data
│   ├── patient_james.py             # James Cooper patient data
│   ├── patient_karen.py             # Karen Foster patient data
│   ├── patient_maya.py              # Maya Rodriguez patient data
│   ├── patient_michael.py           # Michael Torres patient data
│   ├── patients_345.py              # Patients 3-5 data
│   ├── patients_additional.py       # Additional patient data
│   ├── patients_david_enhanced.py   # David Park enhanced data
│   ├── patients_dismissals.py       # Patient dismissal narratives
│   ├── patients_linda_enhanced.py   # Linda Chen enhanced data
│   ├── patients_med_lists.py        # Medication lists
│   ├── patients_rachel_enhanced.py  # Rachel Thompson enhanced data
│   └── patients_referrals.py        # Referral documents
├── tests/
│   ├── conftest.py
│   ├── test_autoimmune.py           # Core agent tests
│   ├── test_collections.py          # Collection schema + manager tests
│   ├── test_diagnostic_engine.py    # Diagnostic engine tests
│   ├── test_rag_engine.py           # RAG engine tests
│   ├── test_export.py               # Export module tests
│   ├── test_api.py                  # API endpoint tests
│   ├── test_timeline_builder.py     # Timeline builder tests
│   └── test_production_readiness.py # Production readiness checks
├── docs/                            # Internal documentation
├── logs/                            # Service logs
├── docker-compose.yml               # autoimmune-api + autoimmune-ui services
├── Dockerfile                       # Multi-stage Python image
├── requirements.txt                 # Python dependencies
├── pyproject.toml                   # Project configuration
├── run.sh                           # Startup script
└── LICENSE                          # Apache 2.0
```

---

## 4. Docker Compose Services

### Service Architecture

```yaml
# docker-compose.yml
services:
  autoimmune-api:
    build: .
    command: uvicorn api.main:app --host 0.0.0.0 --port 8532
    ports: ["8532:8532"]
    environment:
      - AUTO_ANTHROPIC_API_KEY=${AUTO_ANTHROPIC_API_KEY}
      - AUTO_MILVUS_HOST=milvus-standalone
    depends_on: [milvus-standalone]

  autoimmune-ui:
    build: .
    command: streamlit run app/autoimmune_ui.py --server.port 8531
    ports: ["8531:8531"]
    environment:
      - AUTO_ANTHROPIC_API_KEY=${AUTO_ANTHROPIC_API_KEY}
    depends_on: [autoimmune-api]
```

### Startup Sequence

```bash
# 1. Configure
cp .env.example .env
# Edit .env: set AUTO_ANTHROPIC_API_KEY

# 2. Launch
docker compose up -d

# 3. Create collections (first run)
docker compose exec autoimmune-api python3 scripts/setup_collections.py

# 4. Ingest demo data
curl -X POST http://localhost:8532/ingest/demo-data

# 5. Verify
curl http://localhost:8532/health
# → {"status": "healthy"}
```

---

## 5. Milvus Collection Schemas

### 14 Collections Overview

| # | Collection Name | Purpose |
|---|---|---|
| 1 | autoimmune_clinical_documents | Ingested patient clinical records (PDFs) |
| 2 | autoimmune_patient_labs | Lab results with reference range analysis |
| 3 | autoimmune_autoantibody_panels | Autoantibody reference with disease associations |
| 4 | autoimmune_hla_associations | HLA allele to disease risk mapping |
| 5 | autoimmune_disease_criteria | ACR/EULAR classification criteria |
| 6 | autoimmune_disease_activity | Activity scoring systems and thresholds |
| 7 | autoimmune_flare_patterns | Flare prediction biomarker patterns |
| 8 | autoimmune_biologic_therapies | Biologic drug database with PGx |
| 9 | autoimmune_pgx_rules | Pharmacogenomic dosing rules |
| 10 | autoimmune_clinical_trials | Autoimmune clinical trials |
| 11 | autoimmune_literature | Published autoimmune literature |
| 12 | autoimmune_patient_timelines | Diagnostic timeline events |
| 13 | autoimmune_cross_disease | Cross-disease overlap syndromes |
| 14 | genomic_evidence | Shared (read-only) from Stage 2 |

### Unified Schema Pattern

Every collection follows this pattern:

```python
# Primary key
FieldSchema("id", DataType.VARCHAR, is_primary=True, max_length=128)

# Embedding vector (indexed)
FieldSchema("embedding", DataType.FLOAT_VECTOR, dim=384)
# Index: IVF_FLAT, metric: COSINE, nlist: 1024, nprobe: 16

# Shared text field
FieldSchema("text_chunk", DataType.VARCHAR, max_length=3000)

# Collection-specific fields (VARCHAR, INT64, FLOAT)
```

### AutoimmuneCollectionManager

```python
class AutoimmuneCollectionManager:
    """Manages Milvus collections for the Autoimmune Agent."""

    def __init__(self, host="localhost", port=19530, embedding_dim=384):
        # Connects to Milvus with named alias "autoimmune_agent"

    def connect(self) / disconnect(self):
        # Connection lifecycle with retry

    def create_collection(self, name, drop_existing=False) -> Collection:
        # Creates single collection with IVF_FLAT index

    def create_all_collections(self, drop_existing=False) -> Dict[str, Collection]:
        # Creates all 14 collections (skips genomic_evidence if it already exists)

    def insert(self, collection_name, records) -> int:
        # Batch insert with schema-aware field defaults and string truncation

    def insert_batch(self, collection_name, records, batch_size=100) -> int:
        # Batched insertion for large datasets

    def search(self, collection_name, query_embedding, top_k=5, filter_expr=None) -> List[Dict]:
        # Single collection vector search

    def search_all(self, query_embedding, top_k=5, collections=None, max_workers=6) -> Dict[str, List[Dict]]:
        # Parallel search across multiple collections via ThreadPoolExecutor

    def get_collection_stats(self) -> Dict[str, int]:
        # Returns {collection_name: record_count} for all autoimmune collections

    def list_collections(self) -> List[str]:
        # Lists collections matching autoimmune_* or genomic_evidence
```

---

## 6. Pydantic Data Models

### Core Models

```python
class AutoimmuneDisease(str, Enum):
    """13 supported autoimmune diseases."""
    RHEUMATOID_ARTHRITIS = "rheumatoid_arthritis"
    SYSTEMIC_LUPUS = "systemic_lupus_erythematosus"
    MULTIPLE_SCLEROSIS = "multiple_sclerosis"
    TYPE_1_DIABETES = "type_1_diabetes"
    INFLAMMATORY_BOWEL = "inflammatory_bowel_disease"
    PSORIASIS = "psoriasis"
    ANKYLOSING_SPONDYLITIS = "ankylosing_spondylitis"
    SJOGRENS_SYNDROME = "sjogrens_syndrome"
    SYSTEMIC_SCLEROSIS = "systemic_sclerosis"
    MYASTHENIA_GRAVIS = "myasthenia_gravis"
    CELIAC_DISEASE = "celiac_disease"
    GRAVES_DISEASE = "graves_disease"
    HASHIMOTO_THYROIDITIS = "hashimoto_thyroiditis"

class AutoimmunePatientProfile(BaseModel):
    patient_id: str
    age: int  # 0-150
    sex: str  # "M" or "F"
    autoantibody_panel: Optional[AutoantibodyPanel] = None
    hla_profile: Optional[HLAProfile] = None
    biomarkers: Dict[str, float] = {}
    genotypes: Dict[str, str] = {}
    diagnosed_conditions: List[AutoimmuneDisease] = []
    current_medications: List[str] = []
    symptom_duration_months: Optional[int] = None
    family_history: List[str] = []
    # Validator: at least one of autoantibody_panel, hla_profile, or biomarkers required

class AutoantibodyPanel(BaseModel):
    patient_id: str
    collection_date: str  # ISO-8601
    results: List[AutoantibodyResult] = []
    # Properties: positive_antibodies, positive_count

class AutoantibodyResult(BaseModel):
    antibody: str  # "ANA", "anti-dsDNA", "RF", etc.
    value: float
    unit: str = ""
    reference_range: str = ""
    positive: bool = False
    titer: Optional[str] = None  # "1:320"
    pattern: Optional[str] = None  # "homogeneous", "speckled"

class HLAProfile(BaseModel):
    hla_a: List[str] = []
    hla_b: List[str] = []
    hla_c: List[str] = []
    hla_drb1: List[str] = []
    hla_dqb1: List[str] = []
    # Property: all_alleles (flat list of all alleles)

class AutoimmuneAnalysisResult(BaseModel):
    patient_id: str
    disease_activity_scores: List[DiseaseActivityScore] = []
    flare_predictions: List[FlarePredictor] = []
    hla_associations: List[Dict[str, Any]] = []
    biologic_recommendations: List[BiologicTherapy] = []
    critical_alerts: List[str] = []
    cross_agent_findings: List[Dict[str, Any]] = []

class DiseaseActivityScore(BaseModel):
    disease: AutoimmuneDisease
    score_name: str  # "DAS28-CRP", "SLEDAI-2K", etc.
    score_value: float
    level: DiseaseActivityLevel  # REMISSION/LOW/MODERATE/HIGH/VERY_HIGH
    components: Dict[str, float] = {}
    thresholds: Dict[str, float] = {}

class BiologicTherapy(BaseModel):
    drug_name: str
    drug_class: str  # "TNF inhibitor", "IL-6 receptor inhibitor", etc.
    mechanism: str = ""
    indicated_diseases: List[AutoimmuneDisease] = []
    pgx_considerations: List[str] = []
    contraindications: List[str] = []
    monitoring_requirements: List[str] = []
    efficacy_evidence: str = ""
```

### Collection Record Models

Each collection has a record model with `to_embedding_text()`:

```python
class ClinicalDocumentRecord(BaseModel):
    id: str
    text_chunk: str
    patient_id: str = ""
    doc_type: str = ""
    specialty: str = ""
    # ...
    def to_embedding_text(self) -> str:
        parts = [self.text_chunk]
        if self.doc_type: parts.append(f"Document type: {self.doc_type}")
        if self.specialty: parts.append(f"Specialty: {self.specialty}")
        return " ".join(parts)

class LabResultRecord(BaseModel):
    id: str
    test_name: str = ""
    value: float = 0.0
    unit: str = ""
    flag: str = "normal"
    # ...
    def to_embedding_text(self) -> str:
        return f"{self.test_name}: {self.value} {self.unit} ({self.flag}). ..."
```

---

## 7. Configuration Reference

### AutoimmuneSettings

```python
class AutoimmuneSettings(BaseSettings):
    # Environment prefix: AUTO_
    # Loaded from: .env file or environment variables

    model_config = SettingsConfigDict(env_prefix="AUTO_", env_file=".env", extra="ignore")

    # -- Paths --
    PROJECT_ROOT: Path        # Auto-resolved from settings.py location

    # Properties:
    # DATA_DIR -> PROJECT_ROOT / "data"
    # CACHE_DIR -> DATA_DIR / "cache"
    # REFERENCE_DIR -> DATA_DIR / "reference"
    # DEMO_DATA_DIR -> PROJECT_ROOT / "demo_data"

    # -- Milvus --
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # -- 14 Collection Names --
    COLL_CLINICAL_DOCUMENTS: str = "autoimmune_clinical_documents"
    COLL_PATIENT_LABS: str = "autoimmune_patient_labs"
    COLL_AUTOANTIBODY_PANELS: str = "autoimmune_autoantibody_panels"
    COLL_HLA_ASSOCIATIONS: str = "autoimmune_hla_associations"
    COLL_DISEASE_CRITERIA: str = "autoimmune_disease_criteria"
    COLL_DISEASE_ACTIVITY: str = "autoimmune_disease_activity"
    COLL_FLARE_PATTERNS: str = "autoimmune_flare_patterns"
    COLL_BIOLOGIC_THERAPIES: str = "autoimmune_biologic_therapies"
    COLL_PGX_RULES: str = "autoimmune_pgx_rules"
    COLL_CLINICAL_TRIALS: str = "autoimmune_clinical_trials"
    COLL_LITERATURE: str = "autoimmune_literature"
    COLL_PATIENT_TIMELINES: str = "autoimmune_patient_timelines"
    COLL_CROSS_DISEASE: str = "autoimmune_cross_disease"
    COLL_GENOMIC_EVIDENCE: str = "genomic_evidence"  # shared read-only

    # -- Embedding --
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIM: int = 384
    EMBEDDING_BATCH_SIZE: int = 32
    BGE_INSTRUCTION: str = "Represent this sentence for searching relevant passages: "

    # -- LLM --
    ANTHROPIC_API_KEY: str = ""
    LLM_MODEL: str = "claude-sonnet-4-20250514"
    LLM_MAX_TOKENS: int = 4096
    LLM_TEMPERATURE: float = 0.2

    # -- RAG Search --
    TOP_K_PER_COLLECTION: int = 5
    SCORE_THRESHOLD: float = 0.40
    MAX_EVIDENCE_ITEMS: int = 30
    CONVERSATION_MEMORY_SIZE: int = 3

    # -- Collection Weights (sum = 1.00) --
    WEIGHT_CLINICAL_DOCUMENTS: float = 0.18
    WEIGHT_PATIENT_LABS: float = 0.14
    WEIGHT_AUTOANTIBODY_PANELS: float = 0.12
    WEIGHT_HLA_ASSOCIATIONS: float = 0.08
    WEIGHT_DISEASE_CRITERIA: float = 0.08
    WEIGHT_DISEASE_ACTIVITY: float = 0.07
    WEIGHT_FLARE_PATTERNS: float = 0.06
    WEIGHT_BIOLOGIC_THERAPIES: float = 0.06
    WEIGHT_CLINICAL_TRIALS: float = 0.05
    WEIGHT_LITERATURE: float = 0.05
    WEIGHT_PGX_RULES: float = 0.04
    WEIGHT_PATIENT_TIMELINES: float = 0.03
    WEIGHT_CROSS_DISEASE: float = 0.02
    WEIGHT_GENOMIC_EVIDENCE: float = 0.02

    # -- Ports --
    STREAMLIT_PORT: int = 8531
    API_PORT: int = 8532

    # -- Authentication --
    API_KEY: str = ""  # empty = no auth
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8531"
    MAX_REQUEST_SIZE_MB: int = 50  # PDF uploads

    # -- Document Processing --
    MAX_CHUNK_SIZE: int = 2500
    CHUNK_OVERLAP: int = 200
    PDF_DPI: int = 200

    # -- Citation Scoring --
    CITATION_HIGH: float = 0.80
    CITATION_MEDIUM: float = 0.60

    # -- Flare Risk Thresholds --
    FLARE_RISK_IMMINENT: float = 0.8
    FLARE_RISK_HIGH: float = 0.6
    FLARE_RISK_MODERATE: float = 0.4

    # -- Evidence Display --
    MAX_EVIDENCE_TEXT_LENGTH: int = 1500
    MAX_KNOWLEDGE_CONTEXT_ITEMS: int = 25

    # -- Streaming --
    STREAMING_ENABLED: bool = True

    # -- Timeouts --
    REQUEST_TIMEOUT_SECONDS: int = 60
    MILVUS_TIMEOUT_SECONDS: int = 10
    LLM_MAX_RETRIES: int = 3

    # -- Logging --
    LOG_LEVEL: str = "INFO"
    LOG_DIR: str = ""  # empty = PROJECT_ROOT/logs

    # -- Metrics --
    METRICS_ENABLED: bool = True

    # Validator: warns if weights deviate from 1.0 by more than 0.05
    # Property: collection_config -> {name: {weight, label, name}} mapping
    # Property: all_collection_names -> list of all 14 collection names
```

---

## 8. Embedding Strategy

### BGE-small-en-v1.5

| Property | Value |
|---|---|
| Model | BAAI/bge-small-en-v1.5 |
| Dimension | 384 |
| Parameters | 33.4M |
| Max sequence | 512 tokens |
| Metric | Cosine similarity |

### Asymmetric Encoding

```python
# Query embedding (with instruction prefix for asymmetric retrieval)
query_text = "Represent this sentence for searching relevant passages: " + query

# Document embedding (raw text via to_embedding_text(), no prefix)
doc_text = model.to_embedding_text()
```

### to_embedding_text() Pattern

Every record model implements `to_embedding_text()` to produce an optimal embedding string:

```python
# ClinicalDocumentRecord
def to_embedding_text(self) -> str:
    parts = [self.text_chunk]
    if self.doc_type: parts.append(f"Document type: {self.doc_type}")
    if self.specialty: parts.append(f"Specialty: {self.specialty}")
    return " ".join(parts)

# LabResultRecord
def to_embedding_text(self) -> str:
    return f"{self.test_name}: {self.value} {self.unit} ({self.flag}). ..."
```

### Embedding Cache

The RAG engine maintains an in-memory LRU cache (max 256 entries) for query embeddings to avoid redundant encoding of repeated queries.

---

## 9. AutoantibodyInterpreter

### Method Signature

```python
def interpret_autoantibodies(self, panel: AutoantibodyPanel) -> List[Dict[str, Any]]:
```

### 24 Autoantibodies

The knowledge base (`AUTOANTIBODY_DISEASE_MAP`) maps 24 autoantibodies to disease associations:

ANA, anti-dsDNA, anti-Smith, RF, anti-CCP, anti-Scl-70, anti-centromere, anti-SSA/Ro, anti-SSB/La, anti-Jo-1, AChR antibody, anti-tTG IgA, TSI, anti-TPO, anti-RNP, anti-histone, ANCA (c-ANCA/PR3), ANCA (p-ANCA/MPO), anti-Pm-Scl, anti-RNA Polymerase III, anti-cardiolipin, lupus anticoagulant, anti-beta2-glycoprotein I, anti-MuSK.

### Return Structure

```python
[
    {
        "antibody": str,           # e.g., "anti-dsDNA"
        "disease": str,            # e.g., "systemic_lupus_erythematosus"
        "sensitivity": float,      # 0-1
        "specificity": float,      # 0-1
        "value": float,            # measured value
        "titer": str | None,       # e.g., "1:640"
        "pattern": str | None,     # e.g., "homogeneous"
        "note": str,               # clinical note
    }
]
```

### Interpretation Logic

For each positive antibody in the panel, looks up all disease associations from the knowledge base. Returns complete list of findings with sensitivity/specificity for each disease pairing.

---

## 10. HLAAssociationAnalyzer

### Method Signature

```python
def analyze_hla_associations(self, hla_profile: HLAProfile) -> List[Dict[str, Any]]:
```

### 22 HLA Alleles

The knowledge base (`HLA_DISEASE_ASSOCIATIONS`) contains 22 alleles with 30+ disease associations. Notable associations include:

- HLA-B*27:05 -> AS (OR=87.4) -- strongest known HLA-disease association
- HLA-C*06:02 -> Psoriasis (OR=10.0)
- HLA-DQB1*02:01 -> Celiac (OR=7.0)
- HLA-DQB1*03:02 -> T1D (OR=6.5)
- HLA-DRB1*13:01 -> T1D (OR=0.2, PROTECTIVE)

### Matching Logic

1. For each allele in the HLA profile, check exact match (e.g., `HLA-DRB1*04:01`)
2. If no exact match, check broad allele group (e.g., `DRB1*04` matches any `*04:XX` subtype)
3. Sort results by odds ratio (highest risk first)
4. Protective alleles (OR < 1.0) are identified with notes

### Return Structure

```python
[
    {
        "allele": str,           # e.g., "DRB1*04:01"
        "disease": str,          # e.g., "rheumatoid_arthritis"
        "odds_ratio": float,     # e.g., 4.2
        "pmid": str,             # e.g., "20301572"
        "note": str,             # e.g., "Shared epitope hypothesis"
    }
]
```

---

## 11. DiseaseActivityScorer

### Method Signature

```python
def calculate_disease_activity(
    self,
    biomarkers: Dict[str, float],
    conditions: List[AutoimmuneDisease],
) -> List[DiseaseActivityScore]:
```

### 20 Scoring Systems

Defined in `DISEASE_ACTIVITY_THRESHOLDS` knowledge base. Each scoring system includes:
- Disease applicability
- Component list (joints, biomarkers, clinical assessments)
- Threshold cutoffs for remission/low/moderate/high
- Score range
- Reference PMID

### Scoring Logic

1. For each diagnosed condition, find applicable scoring systems
2. Extract CRP and/or ESR from biomarkers
3. Apply simplified marker-based scoring:
   - If CRP/ESR below remission threshold -> REMISSION
   - If below low threshold -> LOW
   - If below moderate threshold -> MODERATE
   - Otherwise -> HIGH
4. Return DiseaseActivityScore with level, components, and thresholds

---

## 12. FlarePredictor

### Method Signature

```python
def predict_flares(
    self,
    biomarkers: Dict[str, float],
    conditions: List[AutoimmuneDisease],
) -> List[FlarePredictor]:
```

### 13 Disease Patterns

Defined in `FLARE_BIOMARKER_PATTERNS`. Each disease has:
- `early_warning_biomarkers`: List of biomarkers to monitor
- `thresholds`: Specific trigger values and windows
- `protective_signals`: Factors that reduce risk

### Risk Calculation

```
base_risk = 0.3
for each early_warning_biomarker present:
    if inflammatory marker (CRP, ESR, IL-6, calprotectin) > 5:  risk += 0.15
    if complement (C3, C4) < 80:                                risk += 0.15
    if albumin < 3.5:                                           risk += 0.10
    else:                                                       mark as protective

risk = clamp(risk, 0.0, 1.0)

IMMINENT if >= 0.8    (configurable via FLARE_RISK_IMMINENT)
HIGH     if >= 0.6    (configurable via FLARE_RISK_HIGH)
MODERATE if >= 0.4    (configurable via FLARE_RISK_MODERATE)
LOW      if < 0.4
```

### Return Structure

```python
FlarePredictor(
    disease=AutoimmuneDisease,
    current_activity=DiseaseActivityLevel,
    predicted_risk=FlareRisk,           # LOW/MODERATE/HIGH/IMMINENT
    risk_score=float,                   # 0.0-1.0
    contributing_factors=List[str],     # e.g., ["Elevated CRP: 12.5"]
    protective_factors=List[str],       # e.g., ["Normal albumin: 4.1"]
    recommended_monitoring=List[str],   # e.g., ["Repeat CRP in 2-4 weeks"]
    time_horizon_days=90,
)
```

---

## 13. BiologicTherapyAdvisor

### Method Signature

```python
def recommend_biologics(
    self,
    conditions: List[AutoimmuneDisease],
    genotypes: Dict[str, str] = None,
) -> List[BiologicTherapy]:
```

### 22 Biologic Therapies

Defined in `BIOLOGIC_THERAPIES` knowledge base. Each therapy includes:
- Drug name, class, and mechanism
- Indicated diseases (list of disease enum values)
- PGx considerations (genotype-specific efficacy/safety data)
- Contraindications
- Monitoring requirements

### Selection Logic

1. For each therapy in the database:
2. Check if any of the patient's diagnosed conditions appear in `indicated_diseases`
3. If match found, include therapy with full PGx, contraindication, and monitoring data
4. Return list of matching therapies

### Drug Classes Covered

TNF inhibitors (5), IL-6R inhibitors (2), Anti-CD20 (2), IL-17A inhibitors (2), IL-12/23 inhibitors (2), IL-23 p19 (1), BLyS inhibitor (1), T-cell modulator (1), JAK inhibitors (3), Integrin inhibitors (2), TYK2 inhibitor (1).

---

## 14. DiagnosticOdysseyAnalyzer

### Classification Criteria Evaluation

```python
def evaluate_classification_criteria(
    self,
    disease: AutoimmuneDisease,
    clinical_data: Dict[str, Any],
) -> Dict[str, Any]:
```

10 criteria sets implemented: RA (2010 ACR/EULAR), SLE (2019 ACR/EULAR), AS (ASAS), SSc (2013 ACR/EULAR), Sjogren's (2016 ACR/EULAR), MS (2017 McDonald), MG (clinical diagnostic), Celiac (ESPGHAN), IBD (Montreal), Psoriasis (clinical).

Returns: total points, threshold, meets_criteria boolean, met/unmet criteria list.

### Overlap Syndrome Detection

```python
def detect_overlap_syndromes(
    self,
    positive_antibodies: List[str],
    diagnosed_conditions: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
```

9 overlap syndromes: MCTD, SLE-RA, Sjogren's-SLE, POTS/hEDS/MCAS triad, SSc-myositis, T1D-celiac, autoimmune thyroid-T1D (APS-2), RA-Sjogren's, lupus-APS.

### Differential Diagnosis

```python
def generate_differential(
    self,
    positive_antibodies: List[str],
    hla_alleles: Optional[List[str]] = None,
    symptoms: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
```

Scoring: antibody specificity x 2.0 + log2(HLA odds ratio) x 0.5. Returns ranked diseases with evidence.

### Diagnostic Odyssey Analysis

```python
def analyze_diagnostic_odyssey(
    self,
    timeline_events: List[Dict[str, Any]],
) -> Dict[str, Any]:
```

Returns: first_symptom_date, diagnosis_date, diagnostic_delay (days/months/years), specialists seen, misdiagnoses, key diagnostic tests.

---

## 15. Multi-Collection RAG Engine

### AutoimmuneRAGEngine

```python
class AutoimmuneRAGEngine:
    def __init__(self, collection_manager, embedder, llm_client, settings=None, knowledge=None):
        # Initializes with all 14 collections and configurable weights
        # Includes: embed cache (256), conversation history (3-turn), disease keyword detection

    def retrieve(self, query, top_k_per_collection=None, ...) -> CrossCollectionResult:
        # Embed query -> parallel search -> weighted merge -> dedup -> knowledge augmentation

    def query(self, question, patient_context=None, **kwargs) -> str:
        # Full RAG: retrieve + synthesize with Claude

    def query_stream(self, question, patient_context=None, **kwargs) -> Generator[str, None, None]:
        # Streaming RAG for real-time UI display

    def search(self, question, **kwargs) -> List[SearchHit]:
        # Evidence-only search without LLM synthesis

    def find_related(self, entity, top_k=5) -> Dict[str, List[SearchHit]]:
        # Cross-collection entity search
```

### Search Pipeline

1. Embed query with BGE-small (with asymmetric query prefix)
2. Detect disease areas from query keywords (12 disease area categories)
3. Build Milvus filter expressions (patient_id, year range) with injection prevention
4. Parallel search across 14 collections via ThreadPoolExecutor (max_workers=6)
5. Apply collection weights and merge results
6. Deduplicate by ID and content hash (MD5 of first 300 chars)
7. Filter by SCORE_THRESHOLD (0.40)
8. Cap at MAX_EVIDENCE_ITEMS (30)
9. Build knowledge context from HLA, autoantibody, therapy, and flare pattern databases
10. Inject conversation history (3-turn)
11. Send to Claude with system prompt and evidence block
12. Return grounded response with citations

### Citation Format

The system prompt instructs Claude to use structured citations:
- `[AutoAb:anti-CCP]` for autoantibodies
- `[HLA:B*27:05]` for HLA alleles
- `[Activity:DAS28-CRP]` for activity scores
- `[Therapy:Adalimumab]` for therapies
- `[Literature:PMID](url)` for published literature
- `[Trial:NCT_ID](url)` for clinical trials

### Collection Weight Table

| Collection | Weight |
|---|---|
| autoimmune_clinical_documents | 0.18 |
| autoimmune_patient_labs | 0.14 |
| autoimmune_autoantibody_panels | 0.12 |
| autoimmune_hla_associations | 0.08 |
| autoimmune_disease_criteria | 0.08 |
| autoimmune_disease_activity | 0.07 |
| autoimmune_flare_patterns | 0.06 |
| autoimmune_biologic_therapies | 0.06 |
| autoimmune_clinical_trials | 0.05 |
| autoimmune_literature | 0.05 |
| autoimmune_pgx_rules | 0.04 |
| autoimmune_patient_timelines | 0.03 |
| autoimmune_cross_disease | 0.02 |
| genomic_evidence | 0.02 |
| **Total** | **1.00** |

### Citation Scoring

| Level | Threshold | Display |
|---|---|---|
| High confidence | >= 0.80 | Full citation |
| Medium confidence | >= 0.60 | Citation with caveat |
| Below threshold | < 0.40 | Filtered out |

---

## 16. Agent Orchestrator

### AutoimmuneAgent

```python
class AutoimmuneAgent:
    def __init__(self, settings=None):
        # Initializes with optional settings

    def analyze_patient(self, profile: AutoimmunePatientProfile) -> AutoimmuneAnalysisResult:
        # Step 1: Interpret autoantibody panel -> disease associations
        # Step 2: Analyze HLA profile -> genetic susceptibility
        # Step 3: Calculate disease activity scores (requires biomarkers + conditions)
        # Step 4: Predict flare risk (requires biomarkers + conditions)
        # Step 5: Recommend biologic therapies with PGx filtering (requires conditions)
        # Step 6: Generate critical alerts from results
        # Return: AutoimmuneAnalysisResult with all findings

    def interpret_autoantibodies(self, panel: AutoantibodyPanel) -> List[Dict]:
        # Maps positive antibodies to AUTOANTIBODY_DISEASE_MAP

    def analyze_hla_associations(self, hla_profile: HLAProfile) -> List[Dict]:
        # Matches alleles against HLA_DISEASE_ASSOCIATIONS

    def calculate_disease_activity(self, biomarkers, conditions) -> List[DiseaseActivityScore]:
        # CRP/ESR-based scoring against DISEASE_ACTIVITY_THRESHOLDS

    def predict_flares(self, biomarkers, conditions) -> List[FlarePredictor]:
        # Pattern matching against FLARE_BIOMARKER_PATTERNS

    def recommend_biologics(self, conditions, genotypes=None) -> List[BiologicTherapy]:
        # Disease-filtered therapy selection from BIOLOGIC_THERAPIES
```

### Alert Generation

Critical alerts are generated for:
- HIGH or VERY_HIGH disease activity scores
- HIGH or IMMINENT flare risk predictions
- HLA associations with odds ratio > 5

---

## 17. Export Pipeline

### Supported Formats

| Format | Function | Output |
|---|---|---|
| Markdown | `to_markdown()` | Structured text with tables |
| PDF | `to_pdf_bytes()` | Styled clinical report (reportlab) |
| FHIR R4 | `to_fhir_r4()` | JSON Bundle |

### FHIR R4 Bundle Structure

```json
{
  "resourceType": "Bundle",
  "type": "collection",
  "timestamp": "ISO-8601",
  "meta": {
    "profile": ["http://hl7.org/fhir/StructureDefinition/Bundle"],
    "source": "HCLS AI Factory - Precision Autoimmune Agent"
  },
  "entry": [
    {"resource": {"resourceType": "Patient", "id": "...", "identifier": [...]}},
    {"resource": {"resourceType": "DiagnosticReport", "status": "final", "code": {...}, "result": [...]}},
    {"resource": {"resourceType": "Observation", "code": {"text": "DAS28-CRP (RA)"}, "valueQuantity": {...}}},
    {"resource": {"resourceType": "Observation", "code": {"text": "Flare Risk Prediction (SLE)"}, ...}}
  ]
}
```

### PDF Export

NVIDIA-branded PDF with:
- Title page with patient ID and timestamp
- Critical alerts section (red text)
- Disease activity scores table (NVIDIA green header)
- Biologic therapy recommendations with PGx
- Clinical query response
- Footer with disclaimer

---

## 18. FastAPI REST Server

### Service Configuration

| Setting | Value |
|---|---|
| Host | 0.0.0.0 |
| Port | 8532 |
| Framework | FastAPI + Uvicorn |
| Auth | Optional API key (AUTO_API_KEY header: X-API-Key) |
| CORS | Configurable origins |
| Request limit | 50 MB (for PDF uploads) |

### Endpoints

| Method | Path | Description |
|---|---|---|
| GET | `/` | Service info |
| GET | `/health` | Detailed health check (Milvus, embedder, LLM, collections, uptime) |
| GET | `/healthz` | Lightweight health probe for landing page |
| POST | `/query` | Full RAG query with Claude synthesis |
| POST | `/query/stream` | Streaming RAG query via SSE |
| POST | `/search` | Evidence-only search (no LLM) |
| POST | `/analyze` | Full autoimmune analysis pipeline |
| POST | `/differential` | Generate differential diagnosis |
| POST | `/ingest/upload` | Upload and ingest a clinical PDF |
| POST | `/ingest/demo-data` | Ingest all demo patient data |
| GET | `/collections` | List collections with stats |
| POST | `/collections/create` | Create all collections |
| POST | `/export` | Export report (markdown, fhir, pdf) |
| GET | `/metrics` | Prometheus-compatible metrics |

### Startup Lifecycle

The API uses FastAPI lifespan with:
1. Centralized logging configuration
2. Milvus connection with retry (2 attempts)
3. BGE-small embedder loading
4. Anthropic LLM client initialization
5. RAG engine, agent, document processor, diagnostic engine, timeline builder initialization
6. Service status banner logging
7. Graceful Milvus disconnect on shutdown

---

## 19. Streamlit UI

### Service Configuration

| Setting | Value |
|---|---|
| Port | 8531 |
| Framework | Streamlit |
| Theme | HCLS AI Factory (NVIDIA green, dark background) |

### 10 Tabs

| Tab | Content |
|---|---|
| 1. Clinical Query | RAG-powered Q&A with evidence citations, collection filtering, streaming |
| 2. Patient Analysis | Full pipeline: antibody + HLA + activity + flare + biologics |
| 3. Document Ingest | Upload clinical PDFs, patient selection, ingestion status |
| 4. Diagnostic Odyssey | Timeline visualization, diagnostic delay, specialists, misdiagnoses |
| 5. Autoantibody Panel | Interactive antibody input, disease association lookup |
| 6. HLA Analysis | HLA allele input, odds ratio display, mechanism notes |
| 7. Disease Activity | Score selection, biomarker input, threshold visualization |
| 8. Flare Prediction | Disease selection, biomarker input, risk scoring |
| 9. Therapy Advisor | Disease-filtered therapy list, PGx details, monitoring |
| 10. Knowledge Base | Collection stats, knowledge version, evidence explorer |

### Cached Service Initialization

Services (collection manager, embedder, LLM client, RAG engine, agent, diagnostic engine, etc.) are initialized once via `@st.cache_resource(ttl=600)` and shared across tabs.

---

## 20. Demo Patient Data

### 9 Patients

| # | Name | Age/Sex | Disease | Key Antibodies | Key HLA | Key Biomarkers |
|---|---|---|---|---|---|---|
| 1 | Sarah Mitchell | 34F | SLE | ANA 1:640 homogeneous, anti-dsDNA, anti-Smith, anti-SSA/Ro | -- | CRP, C3/C4 low, proteinuria |
| 2 | Maya Rodriguez | 28F | POTS/hEDS/MCAS | -- | -- | Tryptase, catecholamines |
| 3 | Linda Chen | 45F | Sjogren's | anti-SSA/Ro, anti-SSB/La, RF | -- | ESR, IgG, complement C4 |
| 4 | David Park | 45M | AS | -- | HLA-B*27 | CRP, ESR |
| 5 | Rachel Thompson | 38F | MCTD | anti-RNP, ANA | -- | CRP, CK |
| 6 | Emma Williams | 34F | MS (RRMS) | -- | HLA-DRB1*15:01 | EDSS, MRI lesion count |
| 7 | James Cooper | 19M | T1D + Celiac | GAD65, IA-2, anti-tTG IgA | HLA-DQ2, HLA-DQ8 | HbA1c, C-peptide, ferritin |
| 8 | Karen Foster | 48F | SSc (dcSSc) | anti-Scl-70 | -- | mRSS, CRP, FVC |
| 9 | Michael Torres | 41M | Graves' Disease | TSI, anti-TPO | -- | TSH, free T4 |

Each patient has a directory under `demo_data/` containing realistic clinical PDF documents (progress notes, lab reports, imaging, pathology, referrals, medication lists) generated by the patient scripts in `scripts/`.

---

## 21. Cross-Agent Integration

### Biomarker Agent Context

```python
def request_biomarker_context(self, patient_id, biomarker_names) -> Dict:
    # Calls Biomarker Agent API for longitudinal biomarker trends
    # Returns: trends for CRP, ESR, complement, etc.
    # Status: stub (returns structured placeholder)
```

### Imaging Agent Context

```python
def request_imaging_context(self, patient_id, body_regions) -> Dict:
    # Calls Imaging Agent API for joint/organ findings
    # Disease-region mapping:
    #   RA -> hands, feet, knees
    #   AS -> spine, sacroiliac joints
    #   SSc -> lungs, heart
    #   SLE -> kidneys, joints
    # Status: stub
```

### Event Publishing

```python
def publish_diagnosis_event(self, patient_id, disease, confidence, evidence) -> Dict:
    # Publishes autoimmune_diagnosis event for other agents
    # Status: stub (production would use Redis/Kafka)
```

### Cross-Agent Enrichment

```python
def enrich_analysis_with_cross_agent(self, result, patient_id) -> AutoimmuneAnalysisResult:
    # Requests biomarker trends for flare contributing factors
    # Requests imaging for diseases with joint/organ involvement
    # Appends findings to result.cross_agent_findings
```

---

## 22. Monitoring & Metrics

### Prometheus Metrics

Available at `GET /metrics` in Prometheus text format.

| Metric | Type | Description |
|---|---|---|
| `autoimmune_agent_up` | Gauge | Whether the agent is running |
| `autoimmune_collection_vectors{collection}` | Gauge | Vector count per collection |
| `autoimmune_agent_uptime_seconds` | Gauge | Agent uptime |

### Health Endpoint

`GET /health` returns:
- `status`: "healthy"
- `milvus_connected`: boolean
- `collections`: count
- `total_vectors`: sum across all collections
- `embedder_loaded`: boolean
- `llm_available`: boolean
- `uptime_seconds`: integer

---

## 23. Testing Strategy

### Test Coverage

| Test File | Coverage |
|---|---|
| test_autoimmune.py | Core agent: analyze_patient, interpret_autoantibodies, HLA analysis, disease activity, flare prediction, biologic recommendations |
| test_collections.py | All 14 collection schemas, field validation, manager methods, search, insert |
| test_diagnostic_engine.py | Classification criteria (10 diseases), overlap syndromes (9), differential diagnosis, odyssey analysis |
| test_rag_engine.py | Engine init, embedding, disease area detection, retrieve, query, streaming, cache |
| test_export.py | Markdown, PDF, FHIR R4 export, bundle structure validation |
| test_api.py | All 14 endpoints, authentication, CORS, error handling |
| test_timeline_builder.py | Timeline construction, event parsing, delay calculation |
| test_production_readiness.py | Settings validation, collection completeness, knowledge integrity |
| **Total** | **455 tests** |

### Running Tests

```bash
python3 -m pytest tests/ -v
# 455 tests passing
```

---

## 24. HCLS AI Factory Integration

### Shared Infrastructure

| Component | Shared With |
|---|---|
| Milvus 2.4 (port 19530) | All agents, RAG pipeline |
| BGE-small-en-v1.5 | All agents |
| Claude API key | All agents |
| `genomic_evidence` collection | Stage 2 RAG pipeline |
| Docker Compose network | All services |

### Agent Family

| Agent | Port (UI/API) | Focus |
|---|---|---|
| CAR-T Intelligence | 8521/8522 | CAR-T cell therapy development |
| Imaging Intelligence | 8523/8524 | Medical imaging analysis |
| Precision Oncology | 8525/8526 | Tumor-specific treatment |
| Precision Biomarker | 8528/8529 | Genomics-informed biomarker interpretation |
| **Precision Autoimmune** | **8531/8532** | **Autoimmune disease intelligence** |

---

## 25. Implementation Sequence

### Quick Start

```bash
# 1. Clone
cd ai_agent_adds/precision_autoimmune_agent

# 2. Install
pip install -r requirements.txt

# 3. Configure
cp .env.example .env
# Set AUTO_ANTHROPIC_API_KEY

# 4. Create all 14 collections
python3 scripts/setup_collections.py

# 5. Run unit tests (455 tests)
python3 -m pytest tests/ -v

# 6. Launch API
uvicorn api.main:app --host 0.0.0.0 --port 8532 &

# 7. Ingest demo patient data
curl -X POST http://localhost:8532/ingest/demo-data

# 8. Launch UI
streamlit run app/autoimmune_ui.py --server.port 8531
```

### Docker Quick Start

```bash
# Edit .env: set AUTO_ANTHROPIC_API_KEY
docker compose up -d

# Create collections (first run)
docker compose exec autoimmune-api python3 scripts/setup_collections.py

# Ingest demo data
curl -X POST http://localhost:8532/ingest/demo-data

# UI: http://localhost:8531
# API: http://localhost:8532/docs
```

---

## Dependencies

```
pydantic>=2.0
pydantic-settings>=2.7
pymilvus>=2.4.0
sentence-transformers>=2.2.0
anthropic>=0.18.0
streamlit>=1.30.0
fastapi>=0.109.0
uvicorn[standard]>=0.27.0
python-dotenv>=1.0.0
loguru>=0.7.0
numpy>=1.24.0
pandas>=2.0.0
plotly>=5.18.0
reportlab>=4.0.0
PyPDF2>=3.0.0
python-multipart>=0.0.6
prometheus-client>=0.20.0
requests>=2.31.0
httpx>=0.25.0
```
