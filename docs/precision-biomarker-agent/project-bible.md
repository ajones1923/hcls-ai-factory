# Precision Biomarker Intelligence Agent — Project Bible

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0
**License:** Apache 2.0

> Complete implementation reference. Import this document as context for Claude Code sessions working on the Precision Biomarker Intelligence Agent.

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
9. [BiologicalAgeCalculator](#9-biologicalagecalculator)
10. [DiseaseTrajectoryAnalyzer](#10-diseasetrajectoryanalyzer)
11. [PharmacogenomicMapper](#11-pharmacogenomicmapper)
12. [CriticalValueEngine](#12-criticalvalueengine)
13. [DiscordanceDetector](#13-discordancedetector)
14. [LabRangeInterpreter](#14-labrangeinterpreter)
15. [GenotypeAdjuster](#15-genotypeadjuster)
16. [Multi-Collection RAG Engine](#16-multi-collection-rag-engine)
17. [Agent Orchestrator](#17-agent-orchestrator)
18. [Export Pipeline](#18-export-pipeline)
19. [FastAPI REST Server](#19-fastapi-rest-server)
20. [Streamlit UI](#20-streamlit-ui)
21. [Sample Patient Data](#21-sample-patient-data)
22. [Monitoring & Metrics](#22-monitoring-metrics)
23. [Testing Strategy](#23-testing-strategy)
24. [HCLS AI Factory Integration](#24-hcls-ai-factory-integration)
25. [Implementation Sequence](#25-implementation-sequence)

---

## 1. Project Overview & Goals

### What This Agent Does

The Precision Biomarker Intelligence Agent transforms raw lab results into clinically actionable, genomics-informed intelligence. Standard lab reports compare values against population-wide reference ranges. This agent goes further by integrating:

- **Patient genotype** (ApoE, MTHFR, CYP variants, HFE, TCF7L2, PNPLA3, DIO2)
- **Pharmacogenomic star alleles** (CYP2D6, CYP2C19, TPMT, DPYD, etc.)
- **Age/sex-stratified reference ranges** (5 brackets, 8 biomarkers)
- **Multi-lab reference comparison** (Quest vs LabCorp vs Function Health optimal)
- **Cross-biomarker discordance detection** (12 pattern rules)
- **Critical value alerting** (21 threshold rules with severity ordering)
- **Epigenetic age estimation** (PhenoAge + GrimAge composite scores)
- **Multi-domain disease risk stratification** (9 disease trajectories)
- **Population-specific screening** (Ashkenazi Jewish carrier panel)

### Key Results

| Metric | Value |
|---|---|
| Total vectors indexed | 652 across 14 Milvus collections |
| Clinical analysis engines | 9 deterministic engines |
| Disease domains | 9 (cardiovascular, diabetes, liver, thyroid, iron, nutritional, kidney, bone, cognitive) |
| PGx genes | 13 (CPIC Level 1A) |
| Critical value rules | 21 (3 severity tiers) |
| Unit tests | 709 passing |
| Demo validation | 65/65 checks passing |
| Export formats | 4 (FHIR R4, PDF, Markdown, CSV) |
| Python files | 57 |
| Lines of code | ~29,000 |

### Pipeline Pattern

```
Patient Biomarkers + Genotypes + Demographics
    │
    ├──> CriticalValueEngine (immediate alerts)
    ├──> BiologicalAgeCalculator (PhenoAge + GrimAge)
    ├──> DiseaseTrajectoryAnalyzer (9 domains)
    ├──> PharmacogenomicMapper (13 genes)
    ├──> GenotypeAdjuster (genotype + age-stratified)
    ├──> DiscordanceDetector (12 patterns)
    ├──> LabRangeInterpreter (3-lab comparison)
    ├──> AJ Carrier Screening (population-specific)
    └──> Age-Stratified Adjustments (5 brackets)
    │
    ▼
RAG Engine (14-collection search + Claude LLM)
    │
    ▼
Export (FHIR R4 | PDF | Markdown | CSV)
```

### HCLS AI Factory Integration

| VAST AI OS Component | Biomarker Agent Role |
|---|---|
| **DataStore** | 16 JSON reference files (652+ records) |
| **DataEngine** | seed_all.py: JSON → BGE-small → Milvus |
| **DataBase** | 14 Milvus collections + 2 sample patients |
| **InsightEngine** | 9 clinical engines + multi-collection RAG |
| **AgentEngine** | PrecisionBiomarkerAgent + Streamlit + FastAPI |

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
- Sufficient for all 9 clinical engines + Milvus + Streamlit simultaneously
- BGE-small inference runs entirely on GPU with <5ms latency

---

## 3. Repository Layout

```
precision_biomarker_agent/
├── src/
│   ├── models.py                    # 20+ Pydantic models (PatientProfile, AnalysisResult, etc.)
│   ├── collections.py               # 14 Milvus collection schemas + BiomarkerCollectionManager
│   ├── rag_engine.py                # BiomarkerRAGEngine (parallel search + Claude)
│   ├── agent.py                     # PrecisionBiomarkerAgent orchestrator
│   ├── biological_age.py            # PhenoAge + GrimAge calculators (~400 lines)
│   ├── disease_trajectory.py        # 9-domain risk analyzer (~1,400 lines)
│   ├── pharmacogenomics.py          # CPIC PGx mapper, 13 genes (~1,400 lines)
│   ├── genotype_adjustment.py       # Genotype + age-stratified adjustments (~1,200 lines)
│   ├── critical_values.py           # Critical value threshold engine (~200 lines)
│   ├── discordance_detector.py      # Cross-biomarker discordance (~300 lines)
│   ├── lab_range_interpreter.py     # Quest vs LabCorp vs optimal (~200 lines)
│   └── export.py                    # FHIR R4 + PDF + Markdown + CSV (~1,000 lines)
├── app/
│   └── biomarker_ui.py              # Streamlit UI, 8 tabs (~1,700 lines)
├── api/
│   └── main.py                      # FastAPI REST server
├── config/
│   └── settings.py                  # PrecisionBiomarkerSettings (Pydantic BaseSettings)
├── data/
│   └── reference/
│       ├── biomarker_reference.json              # 208 biomarker definitions
│       ├── biomarker_genetic_variants.json        # 42 clinically actionable SNPs
│       ├── biomarker_pgx_rules.json              # 29 CPIC pharmacogenomic rules
│       ├── biomarker_disease_trajectories.json   # 39 disease trajectory patterns
│       ├── biomarker_clinical_evidence.json      # 80 clinical evidence records
│       ├── biomarker_nutrition.json              # 50 nutrient-biomarker interactions
│       ├── biomarker_drug_interactions.json      # 51 drug-biomarker effects
│       ├── biomarker_aging_markers.json          # 20 epigenetic clock markers
│       ├── biomarker_genotype_adjustments.json   # 30 genotype-specific adjustments
│       ├── biomarker_monitoring.json             # 30 follow-up protocols
│       ├── biomarker_critical_values.json        # 21 critical value rules
│       ├── biomarker_discordance_rules.json      # 12 discordance patterns
│       ├── biomarker_aj_carrier_screening.json   # 10 AJ carrier screening records
│       ├── biomarker_genomic_evidence.json       # 30 genomic evidence records
│       ├── biomarker_lab_ranges.json             # Multi-lab reference ranges (dict)
│       ├── biomarker_sample_patients.json        # 2 sample patients
│       ├── biomarker_longitudinal_tracking.json  # Longitudinal tracking config (dict)
│       └── nutrition_guidelines_seed.json        # Nutrition guidelines (dict)
├── scripts/
│   ├── seed_all.py                  # Seed all 14 Milvus collections from JSON
│   ├── validate_e2e.py              # End-to-end data layer validation
│   └── demo_validation.py           # 65-check comprehensive demo validation
├── tests/
│   ├── conftest.py                  # Shared fixtures (14 mock collections, 2 patients)
│   ├── test_critical_values.py      # 29 tests
│   ├── test_discordance_detector.py # 25 tests
│   ├── test_lab_range_interpreter.py # 36 tests
│   ├── test_collections.py          # 125 tests
│   ├── test_ui.py                   # 44 tests
│   ├── test_rag_engine.py           # RAG engine tests
│   └── test_edge_cases.py           # Edge case coverage
├── docker-compose.yml               # biomarker-api + biomarker-ui services
├── Dockerfile                       # Multi-stage Python image
├── .env.example                     # Environment variable template
├── requirements.txt                 # Python dependencies
└── LICENSE                          # Apache 2.0
```

---

## 4. Docker Compose Services

### Service Architecture

```yaml
# docker-compose.yml
services:
  biomarker-api:
    build: .
    command: uvicorn api.main:app --host 0.0.0.0 --port 8529
    ports: ["8529:8529"]
    environment:
      - BIOMARKER_ANTHROPIC_API_KEY=${BIOMARKER_ANTHROPIC_API_KEY}
      - BIOMARKER_MILVUS_HOST=milvus-standalone
    depends_on: [milvus-standalone]

  biomarker-ui:
    build: .
    command: streamlit run app/biomarker_ui.py --server.port 8528
    ports: ["8528:8528"]
    environment:
      - BIOMARKER_ANTHROPIC_API_KEY=${BIOMARKER_ANTHROPIC_API_KEY}
    depends_on: [biomarker-api]
```

### Startup Sequence

```bash
# 1. Configure
cp .env.example .env
# Edit .env: set BIOMARKER_ANTHROPIC_API_KEY

# 2. Launch
docker compose up -d

# 3. Seed data (first run only)
docker compose exec biomarker-api python3 scripts/seed_all.py

# 4. Verify
curl http://localhost:8529/health
# → {"status": "healthy"}
```

---

## 5. Milvus Collection Schemas

### 14 Collections Overview

| # | Collection Name | Records | Purpose |
|---|---|---|---|
| 1 | biomarker_reference | 208 | Core biomarker definitions |
| 2 | biomarker_genetic_variants | 42 | Clinically actionable SNPs |
| 3 | biomarker_pgx_rules | 29 | CPIC PGx guidelines |
| 4 | biomarker_disease_trajectories | 39 | Disease risk patterns |
| 5 | biomarker_clinical_evidence | 80 | Published evidence |
| 6 | biomarker_nutrition | 50 | Nutrient interactions |
| 7 | biomarker_drug_interactions | 51 | Medication effects |
| 8 | biomarker_aging_markers | 20 | Epigenetic clock data |
| 9 | biomarker_genotype_adjustments | 30 | Genotype adjustments |
| 10 | biomarker_monitoring | 30 | Follow-up protocols |
| 11 | biomarker_critical_values | 21 | Threshold rules |
| 12 | biomarker_discordance_rules | 12 | Discordance patterns |
| 13 | biomarker_aj_carrier_screening | 10 | AJ carrier screening |
| 14 | genomic_evidence | 30 | Shared (read-only) |

### Unified Schema Pattern

Every collection follows this pattern:

```python
# Primary key
FieldSchema("id", DataType.VARCHAR, is_primary=True, max_length=64)

# Embedding vector (indexed)
FieldSchema("embedding", DataType.FLOAT_VECTOR, dim=384)
# Index: IVF_FLAT, metric: COSINE, nlist: 1024, nprobe: 16

# Text fields (collection-specific)
FieldSchema("name", DataType.VARCHAR, max_length=200)
FieldSchema("description", DataType.VARCHAR, max_length=2000)
# ... additional fields per collection

# Metadata
FieldSchema("category", DataType.VARCHAR, max_length=100)
```

### BiomarkerCollectionManager

```python
class BiomarkerCollectionManager:
    """Manages all 14 Milvus collections for the biomarker agent."""

    def __init__(self, host="localhost", port=19530):
        # Connects to Milvus and initializes all collections

    def create_all(self):
        # Creates all 14 collections with IVF_FLAT indexes

    def insert(self, collection_name, records):
        # Batch insert with BGE-small embedding

    def search(self, collection_name, query_embedding, top_k=5):
        # Single collection search

    def search_all(self, query_embedding, top_k=5):
        # Parallel search across all 14 collections (ThreadPoolExecutor)
```

---

## 6. Pydantic Data Models

### Core Models

```python
class PatientProfile(BaseModel):
    patient_id: str
    age: int
    sex: str  # "M" or "F"
    biomarkers: Dict[str, float]
    genotypes: Optional[Dict[str, str]] = None
    star_alleles: Optional[Dict[str, str]] = None
    medications: Optional[List[str]] = None
    ethnicity: Optional[str] = None

class BiologicalAgeResult(BaseModel):
    chronological_age: int
    biological_age: float
    age_acceleration: float
    phenoage_score: float
    grimage_score: Optional[float] = None

class AnalysisResult(BaseModel):
    patient_profile: PatientProfile
    biological_age: BiologicalAgeResult
    disease_trajectories: List[DiseaseTrajectoryResult] = []
    pgx_results: List[PGxResult] = []
    genotype_adjustments: List[GenotypeAdjustmentResult] = []
    critical_alerts: List[str] = []
    timestamp: str  # ISO-8601
```

### Collection Models

Each collection has a corresponding Pydantic model with a `to_embedding_text()` method:

```python
class BiomarkerReference(BaseModel):
    id: str
    name: str
    category: str
    unit: str
    reference_low: float
    reference_high: float
    clinical_significance: str

    def to_embedding_text(self) -> str:
        return f"{self.name} ({self.category}): {self.clinical_significance}"
```

---

## 7. Configuration Reference

### PrecisionBiomarkerSettings

```python
class PrecisionBiomarkerSettings(BaseSettings):
    # Environment prefix: BIOMARKER_
    # Loaded from: .env file or environment variables

    # ── Paths ──
    PROJECT_ROOT: Path        # Auto-resolved
    DATA_DIR: Path            # PROJECT_ROOT / "data"
    REFERENCE_DIR: Path       # DATA_DIR / "reference"

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # ── 14 Collection Names ──
    COLLECTION_BIOMARKER_REF: str = "biomarker_reference"
    COLLECTION_GENETIC_VARIANTS: str = "biomarker_genetic_variants"
    COLLECTION_PGX_RULES: str = "biomarker_pgx_rules"
    COLLECTION_DISEASE_TRAJECTORIES: str = "biomarker_disease_trajectories"
    COLLECTION_CLINICAL_EVIDENCE: str = "biomarker_clinical_evidence"
    COLLECTION_NUTRITION: str = "biomarker_nutrition"
    COLLECTION_DRUG_INTERACTIONS: str = "biomarker_drug_interactions"
    COLLECTION_AGING_MARKERS: str = "biomarker_aging_markers"
    COLLECTION_GENOTYPE_ADJUSTMENTS: str = "biomarker_genotype_adjustments"
    COLLECTION_MONITORING: str = "biomarker_monitoring"
    COLLECTION_CRITICAL_VALUES: str = "biomarker_critical_values"
    COLLECTION_DISCORDANCE_RULES: str = "biomarker_discordance_rules"
    COLLECTION_AJ_CARRIER_SCREENING: str = "biomarker_aj_carrier_screening"
    COLLECTION_GENOMIC: str = "genomic_evidence"  # Shared (read-only)

    # ── Embeddings ──
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIMENSION: int = 384
    EMBEDDING_BATCH_SIZE: int = 32

    # ── LLM ──
    LLM_PROVIDER: str = "anthropic"
    LLM_MODEL: str = "claude-sonnet-4-6"
    ANTHROPIC_API_KEY: Optional[str] = None

    # ── RAG Search ──
    TOP_K_PER_COLLECTION: int = 5
    SCORE_THRESHOLD: float = 0.4

    # ── Collection Weights (sum = 1.00) ──
    WEIGHT_BIOMARKER_REF: float = 0.12
    WEIGHT_GENETIC_VARIANTS: float = 0.11
    WEIGHT_PGX_RULES: float = 0.10
    WEIGHT_DISEASE_TRAJECTORIES: float = 0.10
    WEIGHT_CLINICAL_EVIDENCE: float = 0.09
    WEIGHT_GENOMIC_EVIDENCE: float = 0.08
    WEIGHT_DRUG_INTERACTIONS: float = 0.07
    WEIGHT_AGING_MARKERS: float = 0.07
    WEIGHT_NUTRITION: float = 0.05
    WEIGHT_GENOTYPE_ADJUSTMENTS: float = 0.05
    WEIGHT_MONITORING: float = 0.05
    WEIGHT_CRITICAL_VALUES: float = 0.04
    WEIGHT_DISCORDANCE_RULES: float = 0.04
    WEIGHT_AJ_CARRIER_SCREENING: float = 0.03

    # ── Services ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8529
    STREAMLIT_PORT: int = 8528

    # ── Citation Scoring ──
    CITATION_HIGH_THRESHOLD: float = 0.75
    CITATION_MEDIUM_THRESHOLD: float = 0.60

    # ── Authentication ──
    API_KEY: str = ""  # Empty = no auth; set for production
```

Weight sum validation: model validator warns if weights deviate from 1.0 by more than 0.05.

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

# Document embedding (raw text, no prefix)
doc_text = model.to_embedding_text()
```

### to_embedding_text() Pattern

Every Pydantic model that maps to a Milvus collection implements `to_embedding_text()`:

```python
# BiomarkerReference
def to_embedding_text(self) -> str:
    return f"{self.name} ({self.category}): {self.clinical_significance}"

# GeneticVariant
def to_embedding_text(self) -> str:
    return f"{self.gene} {self.rsid} {self.risk_allele}: {self.effect_size}"

# PGxRule
def to_embedding_text(self) -> str:
    return f"{self.gene} {self.phenotype} {self.drugs_affected}: {self.recommendation}"
```

---

## 9. BiologicalAgeCalculator

### PhenoAge (Levine 2018)

**9 Input Biomarkers (US clinical units → SI conversion):**

| Biomarker | US Unit | SI Unit | Coefficient |
|---|---|---|---|
| Albumin | g/dL | g/L (×10) | -0.0336 |
| Creatinine | mg/dL | µmol/L (×88.42) | 0.0095 |
| Glucose | mg/dL | mmol/L (÷18.02) | 0.1953 |
| C-Reactive Protein | mg/L | mg/L | 0.0954 |
| Lymphocyte % | % | fraction (÷100) | -0.0120 |
| MCV | fL | fL | 0.0268 |
| RDW | % | % | 0.3306 |
| Alkaline Phosphatase | U/L | U/L | 0.0019 |
| WBC | 10³/µL | 10⁹/L | 0.0554 |

**Calculation:**

```
xb = Σ(coefficient × SI_value) + (0.0804 × chronological_age)
mortality_score = 1 - exp(-exp(xb) × (exp(120 × 0.0076927) - 1) / 0.0076927)
biological_age = 141.50225 + ln(-0.00553 × ln(1 - mortality_score)) / 0.090165
age_acceleration = biological_age - chronological_age
```

**Return structure:**
```python
{
    "phenoage": {
        "chronological_age": float,
        "biological_age": float,
        "age_acceleration": float,
        "mortality_score": float,
        "mortality_risk": str,      # "LOW" | "NORMAL" | "MODERATE" | "HIGH"
        "risk_confidence": str,     # "low" | "moderate" | "high"
        "confidence_interval": {
            "lower": float, "upper": float,
            "confidence_level": 0.95,
            "standard_error": float,
            "note": str
        },
        "top_aging_drivers": [
            {"biomarker": str, "value": float, "si_value": float,
             "coefficient": float, "contribution": float, "direction": str}
        ],
        "missing_biomarkers": [str]
    },
    "grimage": { ... } | None,     # None if no plasma markers
    "biological_age": float,        # shorthand
    "age_acceleration": float,      # shorthand
    "mortality_risk": str           # shorthand
}
```

### GrimAge (Lu 2019)

**6 Plasma Protein Surrogates:**

| Marker | Reference Max | Weight |
|---|---|---|
| GDF-15 | 1,200 pg/mL | 0.25 |
| Cystatin C | 1.0 mg/L | 0.20 |
| Leptin | 10.0 ng/mL | 0.15 |
| PAI-1 | 40 ng/mL | 0.20 |
| TIMP-1 | 250 ng/mL | 0.10 |
| Adrenomedullin | 50 pg/mL | 0.10 |

Returns `None` when no plasma markers are available in the patient data.

---

## 10. DiseaseTrajectoryAnalyzer

### Method Signature

```python
def analyze_all(
    self,
    biomarkers: Dict[str, float],
    genotypes: Dict[str, str],
    age: Optional[float] = None,
    sex: str = "male",
) -> List[Dict[str, Any]]:
```

### 9 Disease Domains

Each domain has an independent `_analyze_*` method that returns a dict with:

```python
{
    "disease": str,                    # e.g., "type2_diabetes"
    "display_name": str,               # e.g., "Type 2 Diabetes"
    "risk_level": str,                 # "CRITICAL" | "HIGH" | "MODERATE" | "LOW"
    "stage": str,                      # Disease-specific stage
    "current_markers": Dict[str, float],
    "genetic_risk_factors": [
        {"gene": str, "genotype": str, "risk_alleles": int, "effect": str}
    ],
    "findings": [str],
    "recommendations": [str],
}
```

Results are sorted by risk severity (CRITICAL first).

### Domain-Specific Genotype Integration

| Domain | Genotype | Effect |
|---|---|---|
| Type 2 Diabetes | TCF7L2 rs7903146 | CT=1 risk allele, TT=2 risk alleles → adjusted glucose thresholds |
| Cardiovascular | APOE | E4 carriers → lower LDL threshold |
| Liver | PNPLA3 rs738409 | GG homozygous → elevated NAFLD risk |
| Thyroid | DIO2 rs225014 | AA/GA → impaired T4→T3 conversion |
| Iron | HFE rs1800562 | GA/AA → hereditary hemochromatosis risk |
| Nutritional | MTHFR rs1801133 | CT/TT → impaired folate metabolism |
| Cognitive | APOE E4 | E3/E4=1.5x risk, E4/E4=3x risk → modifiable risk factor count |

---

## 11. PharmacogenomicMapper

### Method Signature

```python
def map_all(
    self,
    star_alleles: Optional[Dict[str, str]] = None,
    genotypes: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
```

### Return Structure

```python
{
    "gene_results": [
        {
            "gene": str,
            "display_name": str,
            "star_alleles": str,         # e.g., "*1/*4"
            "phenotype": str,            # e.g., "Intermediate Metabolizer"
            "affected_drugs": [
                {"drug": str, "action": str, "recommendation": str, "evidence_level": str}
            ],
            "critical_alerts": [str],
        }
    ],
    "critical_alerts": [str],
    "drugs_to_avoid": [str],
    "drugs_to_adjust": [str],
    "drug_interaction_warnings": [{"gene": str, "drugs": [str], "warning": str}],
    "drug_interactions": [{"drug_a": str, "drug_b": str, "severity": str, "effect": str}],
    "genes_analyzed": int,
    "has_critical_findings": bool,
    "guideline_versions": {"gene": "update_date"},
}
```

### 13 Supported Genes

| Gene | Input | Phenotypes |
|---|---|---|
| CYP2D6 | Star alleles | Poor, Intermediate, Normal, Ultra-rapid |
| CYP2C19 | Star alleles | Poor, Intermediate, Normal, Rapid, Ultra-rapid |
| CYP2C9 | Star alleles | Poor, Intermediate, Normal |
| VKORC1 | Genotype rs9923231 | High, Intermediate, Normal sensitivity |
| SLCO1B1 | Genotype rs4149056 | Poor, Intermediate, Normal function |
| TPMT | Star alleles | Poor, Intermediate, Normal |
| DPYD | Star alleles | Poor, Intermediate, Normal |
| MTHFR | Genotype rs1801133 | Homozygous, Heterozygous, Normal |
| HLA-B*57:01 | Genotype | Positive, Negative |
| G6PD | Genotype | Deficient, Intermediate, Normal |
| HLA-B*58:01 | Genotype | Positive, Negative |
| CYP3A5 | Star alleles | Poor, Intermediate, Normal |
| UGT1A1 | Star alleles | Poor, Intermediate, Normal |

### CPIC Guideline Versions

Each gene entry includes the CPIC guideline publication year, PMID, and last update date for audit trail purposes.

---

## 12. CriticalValueEngine

### Method Signature

```python
def check(self, biomarkers: Dict[str, float]) -> List[CriticalValueAlert]:
```

### Alert Object

```python
class CriticalValueAlert:
    biomarker: str          # Which value triggered
    value: float            # Measured result
    threshold: float        # Exceeded threshold
    direction: str          # "high" or "low"
    severity: str           # "critical" | "urgent" | "warning"
    escalation_target: str  # Routing destination
    clinical_action: str    # Required next step
    cross_checks: List[str] # Related biomarkers to verify
    loinc_code: str         # LOINC identifier
```

### Severity Ordering

```
CRITICAL  →  Immediate notification (e.g., Glucose > 400, Potassium > 6.5)
URGENT    →  Within 4 hours (e.g., Sodium < 125, Hemoglobin < 7.0)
WARNING   →  Next clinical visit (e.g., TSH > 10, Creatinine > 2.0)
```

### Alias Map

The engine maps common biomarker field names to canonical rule names:

```python
"glucose" → "Glucose"
"potassium" → "Potassium"
"calcium" → "Calcium (Total)"
"troponin" → "Troponin I (High-Sensitivity)"
"egfr" → "eGFR (CKD-EPI)"
```

---

## 13. DiscordanceDetector

### Method Signature

```python
def check(self, biomarkers: Dict[str, float]) -> List[DiscordanceFinding]:
```

### 12 Discordance Rules

Examples of clinically significant cross-biomarker discordance patterns:

| Pattern | Condition | Clinical Significance |
|---|---|---|
| LDL normal + ApoB elevated | LDL < 130, ApoB > 90 | Particle number discordance — small dense LDL |
| TSH normal + Free T3 low | TSH 0.5-4.0, Free T3 < 2.3 | Subclinical conversion dysfunction |
| Ferritin low + Hemoglobin normal | Ferritin < 30, Hgb > 12 | Pre-anemic iron depletion |
| HbA1c normal + HOMA-IR elevated | HbA1c < 5.7, HOMA-IR > 2.5 | Insulin resistance preceding glucose dysregulation |

---

## 14. LabRangeInterpreter

### Method Signature

```python
def interpret(
    self,
    biomarkers: Dict[str, float],
    sex: str = "Male",
) -> List[RangeComparison]:
```

### Three-Way Comparison

For each biomarker, compares against:

1. **Quest Diagnostics** reference range
2. **LabCorp** reference range
3. **Function Health** optimal range

Status values: `"normal"`, `"low"`, `"high"`, `"unknown"`

### Sex-Specific Lookup

```python
# Try sex-specific first, then fall back
canonical = self._resolve(biomarker_name)
sex_specific = f"{canonical} ({sex})"
range_data = labs.get(sex_specific) or labs.get(canonical)
```

---

## 15. GenotypeAdjuster

### Genotype Adjustments

```python
def adjust_all(
    self,
    biomarkers: Dict[str, float],
    genotypes: Dict[str, str],
    sex: str = "male",
) -> Dict[str, Any]:
```

### Age-Stratified Adjustments

```python
def apply_age_adjustments(
    self,
    biomarkers: Dict[str, float],
    age: int,
    sex: str = "M",
) -> List[Dict[str, Any]]:
```

### 5 Age Brackets

| Bracket | Ages |
|---|---|
| Pediatric | 0-17 |
| Young Adult | 18-39 |
| Middle Adult | 40-59 |
| Older Adult | 60-79 |
| Elderly | 80+ |

### 8 Age-Stratified Biomarkers

```python
AGE_STRATIFIED_RANGES = {
    "creatinine": {
        "18-39": {"M": (0.7, 1.2), "F": (0.5, 0.9)},
        "40-59": {"M": (0.7, 1.3), "F": (0.6, 1.0)},
        "60-79": {"M": (0.8, 1.4), "F": (0.6, 1.1)},
        "80+":   {"M": (0.9, 1.6), "F": (0.7, 1.3)},
    },
    "egfr": { ... },
    "tsh": { ... },
    "fasting_glucose": { ... },
    "total_cholesterol": { ... },
    "alkaline_phosphatase": { ... },
    "ferritin": { ... },
    "psa": { ... },  # Male only, 40+
}
```

### Return Structure

```python
[
    {
        "biomarker": str,
        "value": float,
        "standard_range": (float, float),
        "age_adjusted_range": (float, float),
        "flag": str,       # "normal" | "low" | "high" | "reclassified"
        "note": str,       # Clinical explanation
    }
]
```

---

## 16. Multi-Collection RAG Engine

### BiomarkerRAGEngine

```python
class BiomarkerRAGEngine:
    def __init__(self, collection_manager, embedder, llm_client, knowledge=None):
        # Initializes with all 14 collections and configurable weights

    def query(self, question: str, patient_profile=None, **kwargs) -> str:
        # Full RAG pipeline: embed → parallel search → augment → Claude generate

    def query_stream(self, question: str, patient_profile=None, **kwargs):
        # Streaming version for real-time UI display
```

### Search Pipeline

1. Embed query with BGE-small (with asymmetric query prefix)
2. Parallel search across 14 collections via ThreadPoolExecutor
3. Weight and merge results by collection weight
4. Filter by SCORE_THRESHOLD (0.4)
5. Augment with knowledge graph context
6. Inject patient profile into prompt
7. Send to Claude Sonnet 4.6
8. Return grounded response with citation scores

### Collection Configuration

Each collection is configured in `config/settings.py` with:
- Collection name (Milvus collection identifier)
- Search weight (0.03 - 0.12, sum = 1.00)
- Top-K per collection (default: 5)

---

## 17. Agent Orchestrator

### PrecisionBiomarkerAgent

```python
class PrecisionBiomarkerAgent:
    def __init__(self):
        self.bio_age = BiologicalAgeCalculator()
        self.trajectory = DiseaseTrajectoryAnalyzer()
        self.pgx = PharmacogenomicMapper()
        self.adjuster = GenotypeAdjuster()
        self.critical = CriticalValueEngine()
        self.discordance = DiscordanceDetector()
        self.lab_ranges = LabRangeInterpreter()
        self.rag_engine = BiomarkerRAGEngine(...)

    def run(self, patient_profile: PatientProfile) -> AgentResponse:
        # Step 1: Critical value check (immediate alerts)
        # Step 2: Biological age calculation
        # Step 3: Disease trajectory analysis (9 domains)
        # Step 4: Pharmacogenomic mapping
        # Step 5: Genotype adjustments
        # Step 6: Discordance detection
        # Step 7: Lab range interpretation
        # Step 8: RAG-augmented synthesis
        # Step 9: Age-stratified adjustments
        # Return: AgentResponse with all results + critical alerts
```

---

## 18. Export Pipeline

### Supported Formats

| Format | Function | Output |
|---|---|---|
| FHIR R4 | `export_fhir_diagnostic_report()` | JSON string (Bundle) |
| PDF | `export_pdf()` | Binary PDF (reportlab) |
| Markdown | `export_markdown()` | Markdown string |
| CSV | `export_csv()` | CSV string |

### FHIR R4 Bundle Structure

```json
{
  "resourceType": "Bundle",
  "type": "diagnostic-report",
  "entry": [
    {"resource": {"resourceType": "Patient", "id": "...", "identifier": [...]}},
    {"resource": {"resourceType": "DiagnosticReport", "status": "final", ...}},
    {"resource": {"resourceType": "Observation", ...}},
    ...
  ]
}
```

### FHIR Validation

```python
def validate_fhir_bundle(bundle: dict) -> List[str]:
    # Checks:
    # 1. Bundle resourceType
    # 2. Entry list presence
    # 3. DiagnosticReport required fields
    # 4. Observation required fields
    # 5. Patient identifier
    # 6. Reference integrity (all refs resolve within bundle)
```

---

## 19. FastAPI REST Server

### Service Configuration

| Setting | Value |
|---|---|
| Host | 0.0.0.0 |
| Port | 8529 |
| Framework | FastAPI + Uvicorn |
| Auth | Optional API key (BIOMARKER_API_KEY) |

### Endpoints

| Method | Path | Description |
|---|---|---|
| GET | `/health` | Health check |
| POST | `/analyze` | Full biomarker analysis |
| POST | `/query` | RAG query with patient context |
| POST | `/export/fhir` | FHIR R4 export |

---

## 20. Streamlit UI

### Service Configuration

| Setting | Value |
|---|---|
| Port | 8528 |
| Framework | Streamlit |
| Theme | HCLS AI Factory (NVIDIA green) |

### 8 Tabs

| Tab | Content |
|---|---|
| 1. Biomarker Analysis | Main analysis with critical alerts, discordance, lab ranges |
| 2. Biological Age | PhenoAge + GrimAge visualization with aging drivers |
| 3. Disease Trajectories | 9-domain risk cards with risk levels and recommendations |
| 4. Pharmacogenomics | Star allele + genotype input, drug recommendations |
| 5. Evidence Explorer | RAG search across 14 collections with source filtering |
| 6. Export | FHIR R4, PDF, Markdown, CSV download |
| 7. Genotype Adjustments | Genotype-aware ranges + age-stratified comparisons |
| 8. Patient Data | Sample patient browser with full clinical context |

### Sample Patient Auto-Load

The `_load_sample_patient()` function maps biomarker JSON keys to Streamlit widget session state keys, pre-filling all input fields for demo flow.

---

## 21. Sample Patient Data

### Patient 1: Male, 45, Ashkenazi Jewish

**Key biomarkers:**

| Biomarker | Value | Unit |
|---|---|---|
| LDL-C | 138 | mg/dL |
| HDL-C | 52 | mg/dL |
| HbA1c | 5.6 | % |
| Fasting Glucose | 98 | mg/dL |
| Fasting Insulin | 8.2 | µIU/mL |
| HOMA-IR | 1.98 | — |
| hs-CRP | 1.8 | mg/L |
| TSH | 2.1 | mIU/L |
| Ferritin | 145 | ng/mL |
| Vitamin D | 38 | ng/mL |
| Lp(a) | 85 | nmol/L |
| Omega-3 Index | 5.8 | % |

**Key genotypes:** APOE E3/E4, MTHFR CT, TCF7L2 CT, PNPLA3 CC, HFE GA, DIO2 GA

**Key star alleles:** CYP2D6 *1/*4, CYP2C19 *1/*2, TPMT *1/*1

### Patient 2: Female, 38, Ashkenazi Jewish

**Key biomarkers:**

| Biomarker | Value | Unit |
|---|---|---|
| LDL-C | 112 | mg/dL |
| HDL-C | 68 | mg/dL |
| HbA1c | 5.2 | % |
| Ferritin | 28 | ng/mL |
| Vitamin D | 32 | ng/mL |
| TSH | 3.4 | mIU/L |
| Omega-3 Index | 4.9 | % |

**Key genotypes:** MTHFR CC, TCF7L2 TT, PNPLA3 CG, DIO2 GG

**Clinical context:** BRCA1 untested (URGENT), active preconception planning (12-18 months), GBA carrier risk 50%

---

## 22. Monitoring & Metrics

### Prometheus Metrics

Enabled via `METRICS_ENABLED=True` in settings.

| Metric | Type | Description |
|---|---|---|
| `biomarker_analysis_duration_seconds` | Histogram | Full analysis latency |
| `biomarker_rag_query_duration_seconds` | Histogram | RAG query latency |
| `biomarker_critical_alerts_total` | Counter | Critical alerts triggered |
| `biomarker_collection_search_latency` | Histogram | Per-collection search latency |

---

## 23. Testing Strategy

### Test Coverage

| Test File | Tests | Coverage |
|---|---|---|
| test_critical_values.py | 29 | All 21 rules, aliases, severity ordering |
| test_discordance_detector.py | 25 | All 12 patterns, edge cases |
| test_lab_range_interpreter.py | 36 | All labs, sex-specific, aliases |
| test_collections.py | 125 | Schemas, models, manager, search, insert |
| test_ui.py | 44 | Source parsing, page config, tabs, helpers |
| test_rag_engine.py | varies | Engine init, weight validation, search |
| test_edge_cases.py | varies | Null handling, empty inputs, boundary values |
| conftest.py | — | 14 mock collections, 2 sample patients |
| **Total** | **709** | |

### Demo Validation Script

```bash
python3 scripts/demo_validation.py
```

65 checks across all 8 tabs:
- Critical values (3 checks)
- Discordance detection (2 checks)
- Lab range interpretation (3 checks)
- Biological age — PhenoAge + GrimAge (9 checks)
- Disease trajectories — 9 domains (13 checks)
- Pharmacogenomics — star alleles + genotypes (8 checks)
- RAG engine initialization (2 checks)
- FHIR export + validation (6 checks)
- Genotype + age adjustments (7 checks)
- Sample patient data integrity (12 checks)

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
| **Precision Biomarker** | **8528/8529** | **Genomics-informed biomarker interpretation** |
| Precision Autoimmune | 8530/8531 | Autoimmune disease intelligence |

---

## 25. Implementation Sequence

### Quick Start

```bash
# 1. Clone
cd ai_agent_adds/precision_biomarker_agent

# 2. Install
pip install -r requirements.txt

# 3. Configure
cp .env.example .env
# Set BIOMARKER_ANTHROPIC_API_KEY

# 4. Seed all 14 collections
python3 scripts/seed_all.py

# 5. Validate data layer
python3 scripts/validate_e2e.py

# 6. Run demo validation (65 checks)
python3 scripts/demo_validation.py

# 7. Run unit tests (709 tests)
python3 -m pytest tests/ -v

# 8. Launch UI
streamlit run app/biomarker_ui.py --server.port 8528

# 9. Launch API (separate terminal)
uvicorn api.main:app --host 0.0.0.0 --port 8529
```

### Docker Quick Start

```bash
# Edit .env: set BIOMARKER_ANTHROPIC_API_KEY
docker compose up -d

# Seed data (first run)
docker compose exec biomarker-api python3 scripts/seed_all.py

# UI: http://localhost:8528
# API: http://localhost:8529/docs
```

---

## Dependencies

```
anthropic>=0.40.0,<1.0
fastapi>=0.104.0,<1.0
uvicorn>=0.24.0,<1.0
streamlit>=1.30.0,<2.0
pydantic>=2.5.0,<3.0
pydantic-settings>=2.1.0,<3.0
pymilvus>=2.3.0,<3.0
sentence-transformers>=2.2.0,<3.0
loguru>=0.7.0,<1.0
reportlab>=4.0.0,<5.0
pandas>=2.0.0,<3.0
plotly>=5.18.0,<6.0
```
