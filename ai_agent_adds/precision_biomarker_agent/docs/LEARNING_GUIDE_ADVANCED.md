# Precision Biomarker Intelligence Agent -- Learning Guide (Advanced)

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0
**Audience:** Engineers extending the agent, adding collections, writing new analysis modules, or deploying to production.

---

## Table of Contents

1. [Prerequisites](#chapter-1-prerequisites)
2. [Deep Dive into the RAG Engine](#chapter-2-deep-dive-into-the-rag-engine)
3. [Vector Search Internals](#chapter-3-vector-search-internals)
4. [Adding a New Collection](#chapter-4-adding-a-new-collection)
5. [The Pharmacogenomics Engine Deep Dive](#chapter-5-the-pharmacogenomics-engine-deep-dive)
6. [Biological Age Algorithms](#chapter-6-biological-age-algorithms)
7. [Disease Trajectory Prediction](#chapter-7-disease-trajectory-prediction)
8. [Genotype-Based Reference Ranges](#chapter-8-genotype-based-reference-ranges)
9. [Clinical Intelligence Modules](#chapter-9-clinical-intelligence-modules)
10. [Export System Deep Dive](#chapter-10-export-system-deep-dive)
11. [Testing Strategies](#chapter-11-testing-strategies)
12. [The Autonomous Agent Pipeline](#chapter-12-the-autonomous-agent-pipeline)
13. [Production Deployment](#chapter-13-production-deployment)
14. [Future Architecture](#chapter-14-future-architecture)

**Appendices:**
- [A. Complete API Reference](#appendix-a-complete-api-reference)
- [B. Configuration Reference](#appendix-b-configuration-reference)
- [C. Collection Schema Reference](#appendix-c-collection-schema-reference)

---

## Chapter 1: Prerequisites

### 1.1 Required Knowledge

Before working with this codebase you should be comfortable with:

- **Python 3.10+** -- Pydantic v2, dataclasses, type hints, async/await.
- **Vector databases** -- Milvus, approximate nearest-neighbor search, IVF indices.
- **Embeddings** -- Sentence Transformers, BAAI/bge-small-en-v1.5, cosine similarity.
- **Clinical genomics** -- Star alleles, CPIC guidelines, pharmacogenomics, VCF format.
- **FastAPI** -- Dependency injection, lifespan events, middleware, Pydantic schemas.
- **Docker** -- Multi-stage builds, compose networking, health checks.

### 1.2 Codebase Map

The agent lives at `ai_agent_adds/precision_biomarker_agent/` within the HCLS AI Factory monorepo. Every source file and its line count is listed below.

#### Source Modules (`src/`)

| File                      | Lines | Purpose                                           |
|---------------------------|------:|---------------------------------------------------|
| `pharmacogenomics.py`     | 1,503 | Star allele to metabolizer phenotype mapping (CPIC)|
| `disease_trajectory.py`   | 1,421 | Pre-symptomatic disease trajectory prediction      |
| `collections.py`          | 1,391 | Milvus collection schemas and manager              |
| `knowledge.py`            | 1,326 | Static knowledge graph (domains, PGx, PhenoAge)    |
| `export.py`               | 1,392 | Markdown, JSON, PDF, CSV, FHIR R4 export           |
| `genotype_adjustment.py`  | 1,225 | Genotype- and age-stratified reference ranges       |
| `report_generator.py`     |   993 | 12-section clinical report generation               |
| `models.py`               |   786 | Pydantic models and enums for all data structures   |
| `agent.py`                |   610 | Autonomous agent pipeline (plan/analyze/search)     |
| `rag_engine.py`           |   573 | Multi-collection RAG engine                         |
| `discordance_detector.py` |   299 | Cross-biomarker discordance detection               |
| `lab_range_interpreter.py`|   221 | Standard vs optimal range interpretation            |
| `translation.py`          |   217 | Multi-language report translation                   |
| `critical_values.py`      |   179 | Critical value threshold checking                   |
| `audit.py`                |    83 | Audit logging for PHI access                        |
| `__init__.py`             |     1 | Package marker                                      |
| **Total**                 |**12,628**|                                                  |

#### API Layer (`api/`)

| File                    | Lines | Purpose                                    |
|-------------------------|------:|--------------------------------------------|
| `main.py`               |   465 | FastAPI app, lifespan, middleware, core endpoints |
| `routes/analysis.py`    |   ~300| `/v1/analyze`, `/v1/biological-age`, `/v1/pgx`, `/v1/query` |
| `routes/reports.py`     |   ~250| `/v1/report/generate`, `/v1/report/{id}/pdf`, FHIR export   |
| `routes/events.py`      |   ~200| Cross-modal event ingestion and alert dispatch              |

#### Application Layer (`app/`)

| File                | Lines | Purpose                                     |
|---------------------|------:|---------------------------------------------|
| `biomarker_ui.py`   | 1,863 | Streamlit UI (port 8528)                    |
| `patient_360.py`    |   670 | Patient 360-degree dashboard                |
| `protein_viewer.py` |   168 | 3D protein structure viewer                 |

#### Configuration (`config/`)

| File           | Lines | Purpose                                        |
|----------------|------:|------------------------------------------------|
| `settings.py`  |   139 | Pydantic BaseSettings with `BIOMARKER_` prefix |

#### Tests (`tests/`)

| File                           | Lines | Test count |
|--------------------------------|------:|-----------:|
| `test_edge_cases.py`           |   972 |         69 |
| `test_api.py`                  | 1,080 |         59 |
| `test_disease_trajectory.py`   |   509 |         48 |
| `test_export.py`               |   453 |         46 |
| `test_ui.py`                   |   610 |         39 |
| `test_models.py`               |   585 |         39 |
| `test_lab_range_interpreter.py`|   460 |         37 |
| `test_biological_age.py`       |   406 |         30 |
| `test_critical_values.py`      |   390 |         28 |
| `test_pharmacogenomics.py`     |   380 |         27 |
| `test_genotype_adjustment.py`  |   332 |         26 |
| `test_discordance_detector.py` |   378 |         25 |
| `test_collections.py`          |   279 |         22 |
| `test_report_generator.py`     |   348 |         21 |
| `test_rag_engine.py`           |   273 |         21 |
| `test_integration.py`          |   540 |         21 |
| `test_longitudinal.py`         |   162 |         18 |
| `test_agent.py`                |   307 |         16 |
| `conftest.py`                  |   307 |          - |
| **Total**                      |**8,772**| **709** |

### 1.3 Key Dependencies

```
sentence-transformers   # BAAI/bge-small-en-v1.5, 384-dim embeddings
pymilvus                # Milvus Python SDK
anthropic               # Claude API client
fastapi / uvicorn       # REST API server
streamlit               # Interactive UI
pydantic / pydantic-settings  # Configuration and data models
reportlab               # PDF report generation (Platypus engine)
loguru                   # Structured logging
```

### 1.4 Port Assignments

| Service         | Port |
|-----------------|------|
| Streamlit UI    | 8528 |
| FastAPI API     | 8529 |
| Milvus          | 19530|

---

## Chapter 2: Deep Dive into the RAG Engine

**File:** `src/rag_engine.py` (573 lines)

### 2.1 Architecture Overview

The `BiomarkerRAGEngine` class implements a multi-collection Retrieval-Augmented Generation pipeline. It searches across all 14 Milvus collections simultaneously using a `ThreadPoolExecutor` (delegated to the collection manager), merges results with knowledge graph context, and generates grounded LLM responses via Claude.

```
User Question
    |
    v
[1] Embed query (BGE-small-en-v1.5, 384 dims)
    |
    v
[2] Determine collections to search (14 total, or filtered subset)
    |
    v
[3] Build per-collection filter expressions
    |   - Disease area filter (diabetes, cardiovascular, liver, ...)
    |   - Year range filter (clinical evidence only)
    |
    v
[4] Parallel search across all collections (ThreadPoolExecutor)
    |
    v
[5] Deduplicate + Citation scoring + Rank by weighted score
    |
    v
[6] Knowledge graph augmentation (domains, PGx, PhenoAge, biomarkers)
    |
    v
CrossCollectionResult (max 30 merged hits)
    |
    v
[7] Build prompt with evidence, knowledge context, patient profile
    |
    v
[8] LLM generation (Claude, max_tokens=2048, temperature=0.7)
```

### 2.2 The `retrieve()` Method

This is the core retrieval method. It accepts an `AgentQuery` and returns a `CrossCollectionResult`:

```python
def retrieve(self, query: AgentQuery,
             top_k_per_collection: int = None,
             collections_filter: List[str] = None,
             year_min: int = None,
             year_max: int = None,
             conversation_context: str = None) -> CrossCollectionResult:
```

**Key parameters:**

- `top_k_per_collection`: Max results per collection. Default: `settings.TOP_K_PER_COLLECTION` (5).
- `collections_filter`: Optional list of collection names. If `None`, searches all 14.
- `year_min` / `year_max`: Applied only to `biomarker_clinical_evidence` via the `year` field.
- `conversation_context`: For multi-turn queries; limited to 2,000 chars, prepended to search text.

**Step-by-step flow:**

1. **Embed query** -- Calls `_embed_query()`, which prepends the BGE instruction prefix `"Represent this sentence for searching relevant passages: "` to the question text, then calls `embedder.embed_text()`.
2. **Build filters** -- For collections with `has_disease_area: True`, detects disease area keywords in the question using `_detect_disease_area()`. Filter expressions use Milvus boolean syntax (e.g., `disease_area == "cardiovascular"`). Input is validated with a safe-character regex to prevent injection.
3. **Parallel search** -- Delegates to `collections.search_all()` which uses `ThreadPoolExecutor`. Each collection is searched independently with its own filter expression.
4. **Merge and rank** -- Deduplicates by ID and text prefix (first 200 chars), sorts by weighted score descending, caps at `MAX_MERGED_RESULTS = 30`.

### 2.3 Score Weighting Math

Every search hit receives a weighted score that combines the raw cosine similarity with the collection's importance weight:

```python
weighted_score = min(raw_score * (1 + weight), 1.0)
```

Where `weight` is the collection-specific weight from settings. This formula provides a bounded boost: a collection with weight 0.12 boosts scores by up to 12%. The `min(..., 1.0)` clamp prevents scores from exceeding 1.0.

**Collection weights (must sum to ~1.0):**

| Collection              | Weight | Label           |
|-------------------------|-------:|-----------------|
| `biomarker_reference`   |   0.12 | BiomarkerRef    |
| `genetic_variants`      |   0.11 | GeneticVariant  |
| `pgx_rules`             |   0.10 | PGxRule         |
| `disease_trajectories`  |   0.10 | DiseaseTrajectory|
| `clinical_evidence`     |   0.09 | ClinicalEvidence|
| `genomic_evidence`      |   0.08 | Genomic         |
| `drug_interactions`     |   0.07 | DrugInteraction |
| `aging_markers`         |   0.07 | AgingMarker     |
| `nutrition`             |   0.05 | Nutrition       |
| `genotype_adjustments`  |   0.05 | GenotypeAdj     |
| `monitoring`            |   0.05 | Monitoring      |
| `critical_values`       |   0.04 | CriticalValue   |
| `discordance_rules`     |   0.04 | DiscordanceRule |
| `aj_carrier_screening`  |   0.03 | AJCarrierScreen |
| **Sum**                 |**1.00**|                 |

### 2.4 Citation Relevance Scoring

Each hit is tagged with a relevance level based on the raw similarity score before weighting:

```python
if raw_score >= settings.CITATION_HIGH_THRESHOLD:    # 0.75
    relevance = "high"
elif raw_score >= settings.CITATION_MEDIUM_THRESHOLD:  # 0.60
    relevance = "medium"
else:
    relevance = "low"
```

The relevance tag is injected into the LLM prompt as `[high relevance]`, `[medium relevance]`, or `[low relevance]` next to each citation. The system prompt instructs the LLM to "prioritize [high relevance] citations."

### 2.5 The System Prompt

The system prompt (`BIOMARKER_SYSTEM_PROMPT`) is a 40-line instruction set that defines the agent's nine expertise domains:

1. Biological Aging (PhenoAge, GrimAge, epigenetic clocks)
2. Pre-Symptomatic Disease Detection (trajectories, timelines)
3. Pharmacogenomic Drug-Gene Interactions (CPIC, star alleles)
4. Genotype-Adjusted Reference Ranges (MTHFR, APOE, PNPLA3, etc.)
5. Nutritional Genomics (MTHFR/methylfolate, FADS1/omega-3, VDR/vitamin D)
6. Cardiovascular Risk Stratification (Lp(a), ApoB, APOE, PCSK9)
7. Liver Health Assessment (PNPLA3 I148M, TM6SF2, FIB-4)
8. Iron Metabolism (HFE C282Y/H63D, hemochromatosis)
9. Ashkenazi Jewish Carrier Screening (10-gene AJ panel)

It instructs the LLM to cite evidence using collection labels, specify units, provide genotype-specific interpretation, highlight critical findings, and flag cross-modal triggers.

### 2.6 Prompt Construction

The `_build_prompt()` method assembles the final prompt from four sections:

```
## Retrieved Evidence

### Evidence from BiomarkerRef
1. [BiomarkerRef:albumin] [high relevance] (score=0.892) ...

### Evidence from ClinicalEvidence
1. [ClinicalEvidence:PMID 29676998](https://pubmed.ncbi.nlm.nih.gov/29676998/) ...

### Knowledge Graph Context
PhenoAge Clock Context: ...

### Patient Profile Context
Age: 45, Sex: M
Biomarkers: albumin: 4.1, creatinine: 0.9, ...
Genotypes: rs1801133: CT, ...
Star Alleles: CYP2D6: *1/*4, ...

---

## Question

What does my HbA1c of 5.8% mean given my TCF7L2 CT genotype?

Please provide a comprehensive answer grounded in the evidence above. ...
```

Clinical evidence citations include clickable PubMed URLs: `[ClinicalEvidence:PMID 29676998](https://pubmed.ncbi.nlm.nih.gov/29676998/)`.

### 2.7 Cross-Collection Entity Linking

The `find_related()` method enables cross-collection entity discovery:

```python
engine.find_related("MTHFR")
# Returns: {
#   "biomarker_genetic_variants": [SearchHit(...)],
#   "biomarker_nutrition": [SearchHit(...)],
#   "biomarker_genotype_adjustments": [SearchHit(...)],
# }
```

This powers queries like "show me everything about MTHFR" or "find all CYP2D6 drug interactions" spanning all 14 collections.

---

## Chapter 3: Vector Search Internals

### 3.1 Index Type: IVF_FLAT

All 14 collections use `IVF_FLAT` (Inverted File with Flat quantization) as the index type. This partitions the vector space into clusters using k-means, then performs exhaustive search within the selected clusters.

```python
index_params = {
    "index_type": "IVF_FLAT",
    "metric_type": "COSINE",
    "params": {"nlist": 128}
}
```

- **nlist=128**: Number of Voronoi cells (clusters). At ingest time, each vector is assigned to the nearest of 128 centroids.
- **nprobe=16**: At query time, the 16 nearest clusters are searched. Higher nprobe means better recall at the cost of latency.

The recall/latency tradeoff: with `nprobe=16` out of `nlist=128`, roughly 12.5% of the index is scanned. For biomarker collections (hundreds to low thousands of records), this provides near-perfect recall with sub-millisecond search times.

### 3.2 Distance Metrics: COSINE vs L2 vs IP

The agent uses `COSINE` similarity as its distance metric:

| Metric  | Formula                              | Range     | Use Case                           |
|---------|--------------------------------------|-----------|------------------------------------|
| COSINE  | `1 - cos(A, B)`                     | [0, 2]    | Normalized embeddings (BGE)        |
| L2      | `||A - B||_2`                        | [0, inf)  | Raw distance, sensitive to magnitude|
| IP      | `A . B`                              | (-inf, inf)| Maximizes dot product              |

**Why COSINE?** BGE-small-en-v1.5 produces L2-normalized embeddings, so COSINE and IP are mathematically equivalent. COSINE is chosen because Milvus returns similarity scores in [0, 1] for COSINE, which maps naturally to the citation relevance thresholds (0.75 high, 0.60 medium).

### 3.3 BGE Embedding Model

The agent uses `BAAI/bge-small-en-v1.5`:

- **Dimensions:** 384
- **Model size:** ~33M parameters (~130MB)
- **Sequence length:** 512 tokens max
- **Instruction-tuned:** Uses the prefix `"Represent this sentence for searching relevant passages: "` for queries (but not for documents).

```python
# Query embedding (with instruction prefix)
prefix = "Represent this sentence for searching relevant passages: "
query_vec = embedder.embed_text(prefix + "What affects CYP2D6 metabolism?")

# Document embedding (no prefix)
doc_vec = embedder.embed_text("CYP2D6 is a cytochrome P450 enzyme that...")
```

### 3.4 Search Parameters

```python
search_params = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16}
}
```

The `SCORE_THRESHOLD` setting (default 0.4) filters out hits below minimum relevance. Any hit with `score < 0.4` is discarded before ranking. This prevents low-quality noise from reaching the LLM prompt.

### 3.5 Embedding Pipeline

```
Input Text
    |
    v
SentenceTransformer("BAAI/bge-small-en-v1.5")
    |
    v
model.encode(text)  -->  numpy array (384,)
    |
    v
.tolist()  -->  List[float] (384 elements)
    |
    v
Milvus insert / search
```

At API startup, the model is loaded once and shared across requests:

```python
class _Embedder:
    def __init__(self):
        self.model = SentenceTransformer(settings.EMBEDDING_MODEL)

    def embed_text(self, text: str) -> List[float]:
        return self.model.encode(text).tolist()
```

---

## Chapter 4: Adding a New Collection

This chapter walks through adding a hypothetical `biomarker_microbiome` collection in 10 steps.

### Step 1: Define the Pydantic Model

Add to `src/models.py`:

```python
class MicrobiomeMarker(BaseModel):
    """Microbiome-biomarker interaction -- maps to biomarker_microbiome collection."""
    id: str = Field(..., max_length=100, description="Unique marker identifier")
    organism: str = Field("", max_length=100, description="Bacterial species/genus")
    biomarker_affected: str = Field("", max_length=100, description="Biomarker name")
    mechanism: str = Field("", max_length=2000, description="Mechanism of action")
    text_chunk: str = Field(..., max_length=3000, description="Text for embedding")
    disease_area: str = Field("", max_length=50, description="Disease area tag")
```

### Step 2: Define the Milvus Schema

Add to `src/collections.py`:

```python
MICROBIOME_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="organism", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="biomarker_affected", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="mechanism", dtype=DataType.VARCHAR, max_length=2000),
    FieldSchema(name="text_chunk", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="disease_area", dtype=DataType.VARCHAR, max_length=50),
]
MICROBIOME_SCHEMA = CollectionSchema(
    fields=MICROBIOME_FIELDS,
    description="Microbiome-biomarker interactions",
)
```

### Step 3: Register in BiomarkerCollectionManager

In `collections.py`, add the collection to the schema registry dict (follow the existing pattern for `_COLLECTION_SCHEMAS`):

```python
"biomarker_microbiome": MICROBIOME_SCHEMA,
```

Add the collection name to `__init__` where collections are listed, and add it to `ensure_collections()`.

### Step 4: Add the Weight Setting

In `config/settings.py`:

```python
WEIGHT_MICROBIOME: float = 0.04
```

Adjust other weights so the total still sums to ~1.0. Run the `_validate_settings` model validator to confirm.

### Step 5: Register in COLLECTION_CONFIG

In `src/rag_engine.py`:

```python
"biomarker_microbiome": {
    "weight": settings.WEIGHT_MICROBIOME,
    "label": "Microbiome",
    "has_disease_area": True,
    "year_field": None,
},
```

### Step 6: Add the Setting to env_prefix

The env var is automatically named `BIOMARKER_WEIGHT_MICROBIOME` thanks to the `env_prefix="BIOMARKER_"` in `PrecisionBiomarkerSettings.model_config`.

### Step 7: Create a Seed Script

Create `scripts/seed_microbiome.py` that reads source data, embeds text chunks, and inserts into Milvus:

```python
from sentence_transformers import SentenceTransformer
from pymilvus import Collection

model = SentenceTransformer("BAAI/bge-small-en-v1.5")

records = load_microbiome_data()  # Your data loading function
embeddings = model.encode([r["text_chunk"] for r in records])

collection = Collection("biomarker_microbiome")
collection.insert([
    [r["id"] for r in records],
    embeddings.tolist(),
    [r["organism"] for r in records],
    # ... remaining fields
])
collection.flush()
```

### Step 8: Update conftest.py

Add the new collection to the mock collection manager's `collection_names` list:

```python
collection_names = [
    # ... existing 14 collections ...
    "biomarker_microbiome",
]
```

### Step 9: Write Tests

Create `tests/test_microbiome.py` following the existing test patterns. Test at minimum:
- Schema creation
- Insert and search round-trip
- Weight application in RAG engine
- Disease area filtering

### Step 10: Verify End-to-End

```bash
# Start Milvus
docker compose up -d milvus-standalone

# Seed the new collection
python scripts/seed_microbiome.py

# Run tests
pytest tests/test_microbiome.py -v

# Verify via API
curl http://localhost:8529/collections | python -m json.tool
```

---

## Chapter 5: The Pharmacogenomics Engine Deep Dive

**File:** `src/pharmacogenomics.py` (1,503 lines)

### 5.1 Architecture

The `PharmacogenomicMapper` class implements a pure-computation engine that maps star allele diplotypes to metabolizer phenotypes and drug-specific dosing recommendations. It requires no LLM calls or database queries -- all knowledge is embedded in the `PGX_GENE_CONFIGS` dictionary.

### 5.2 The Fourteen Pharmacogenes

| Gene    | Role                              | CPIC Level | Key Drugs                              |
|---------|-----------------------------------|:----------:|----------------------------------------|
| CYP2D6  | Metabolizes ~25% of drugs         | 1A         | Codeine, tramadol, tamoxifen           |
| CYP2C19 | Proton pump inhibitors, antiplatelets| 1A      | Clopidogrel, omeprazole, voriconazole  |
| CYP2C9  | NSAIDs, warfarin metabolism       | 1A         | Warfarin, celecoxib, phenytoin         |
| CYP3A5  | Immunosuppressant metabolism      | 1A         | Tacrolimus                             |
| SLCO1B1 | Hepatic drug transporter          | 1A         | Simvastatin, atorvastatin              |
| VKORC1  | Warfarin target sensitivity       | 1A         | Warfarin                               |
| MTHFR   | Folate metabolism enzyme          | Info       | Methotrexate (adjunctive)              |
| TPMT    | Thiopurine metabolism             | 1A         | Azathioprine, 6-mercaptopurine         |
| DPYD    | Fluoropyrimidine metabolism       | 1A         | 5-FU, capecitabine                     |

### 5.3 Star Allele to Phenotype Mapping

Each gene has an `allele_to_phenotype` dictionary that maps diplotype strings to phenotype labels:

```python
PGX_GENE_CONFIGS = {
    "CYP2D6": {
        "display_name": "CYP2D6",
        "description": "Cytochrome P450 2D6 -- metabolizes ~25% of drugs",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*4": "Intermediate Metabolizer",
            "*4/*4": "Poor Metabolizer",
            "*1/*1xN": "Ultra-rapid Metabolizer",
            # ... 16 diplotype combinations
        },
        "drug_recommendations": { ... },
    },
    # ... 13 more genes
}
```

### 5.4 Metabolizer Phenotype Classification

The `MetabolizerPhenotype` enum defines four standard CPIC categories:

| Phenotype         | Enum Value      | Clinical Meaning                              |
|-------------------|-----------------|-----------------------------------------------|
| Ultra-rapid       | `ultra_rapid`   | Excess enzyme activity; rapid drug clearance   |
| Normal            | `normal`        | Standard enzyme activity; use standard dosing  |
| Intermediate      | `intermediate`  | Reduced activity; consider dose adjustment     |
| Poor              | `poor`          | Minimal/no activity; avoid or reduce dose      |

Non-CYP genes use specialized terminology:
- **SLCO1B1**: Normal Function / Intermediate Function / Poor Function (transporter activity)
- **VKORC1**: Normal Sensitivity / Intermediate Sensitivity / High Sensitivity (drug target)
- **MTHFR**: Normal Activity / Intermediate Activity / Reduced Activity (enzyme activity)

### 5.5 Drug-Specific Dosing Recommendations

Each drug entry maps every possible phenotype to a structured recommendation:

```python
"codeine": {
    "Poor Metabolizer": {
        "recommendation": "AVOID codeine -- no conversion to morphine, will be ineffective.",
        "action": "AVOID",
        "alert_level": "CRITICAL",
    },
    "Ultra-rapid Metabolizer": {
        "recommendation": "AVOID codeine -- excess conversion to morphine, "
                          "risk of fatal respiratory depression.",
        "action": "AVOID",
        "alert_level": "CRITICAL",
    },
}
```

**Action categories:**
- `STANDARD_DOSING` -- No change needed
- `DOSE_REDUCTION` -- Reduce dose per recommendation
- `DOSE_ADJUSTMENT` -- Adjust dose (up or down)
- `CONSIDER_ALTERNATIVE` -- Current drug may work but alternative preferred
- `AVOID` -- Do not use this drug
- `CONTRAINDICATED` -- Absolute contraindication (FDA/EMA mandated)

**Alert levels:** `INFO`, `WARNING`, `CRITICAL`

### 5.6 CPIC Level Evidence

Every gene entry includes version tracking:

```python
CPIC_GUIDELINE_VERSIONS = {
    "CYP2D6": {"version": "2019", "pmid": "33387367", "update": "2020-12", "level": "1A"},
    "CYP2C19": {"version": "2022", "pmid": "34697867", "update": "2022-12", "level": "1A"},
    # ...
}
```

### 5.7 The `map_all()` Method

```python
pgx_mapper = PharmacogenomicMapper()
results = pgx_mapper.map_all(
    star_alleles={"CYP2D6": "*4/*4", "CYP2C19": "*1/*2"},
    genotypes={"rs1801133": "CT"},  # MTHFR
)
# Returns: {
#   "gene_results": [
#     {"gene": "CYP2D6", "star_alleles": "*4/*4",
#      "phenotype": "Poor Metabolizer", "affected_drugs": [...]},
#     {"gene": "CYP2C19", "star_alleles": "*1/*2",
#      "phenotype": "Intermediate Metabolizer", "affected_drugs": [...]},
#     {"gene": "MTHFR", "genotype": "CT",
#      "phenotype": "Intermediate Activity", "affected_drugs": [...]},
#   ]
# }
```

### 5.8 Adding a New Gene

To add a new pharmacogene (e.g., `NAT2`):

1. Add CPIC version info to `CPIC_GUIDELINE_VERSIONS`.
2. Add the full gene config to `PGX_GENE_CONFIGS` with `allele_to_phenotype` and `drug_recommendations`.
3. Add test cases to `tests/test_pharmacogenomics.py`.
4. The gene is automatically picked up by `map_all()`.

---

## Chapter 6: Biological Age Algorithms

**File:** `src/biological_age.py` (408 lines)

### 6.1 PhenoAge (Levine 2018)

PhenoAge estimates biological age from 9 routine blood biomarkers using a Gompertz mortality model trained on NHANES III data.

**Reference:** Levine et al., "An epigenetic biomarker of aging for lifespan and healthspan", *Aging* 2018; 10(4):573-591. PMID: 29676998.

### 6.2 The Nine Biomarkers and Coefficients

| Biomarker               | Coefficient | Direction   | Units (input) | Units (SI)  |
|--------------------------|------------|-------------|---------------|-------------|
| Albumin                  | -0.0336    | Protective  | g/dL          | g/L         |
| Creatinine               |  0.0095    | Aging       | mg/dL         | umol/L      |
| Glucose                  |  0.1953    | Aging       | mg/dL         | mmol/L      |
| ln(CRP)                  |  0.0954    | Aging       | mg/L (ln)     | ln(mg/L)    |
| Lymphocyte %             | -0.0120    | Protective  | %             | %           |
| MCV                      |  0.0268    | Aging       | fL            | fL          |
| RDW                      |  0.3306    | Aging       | %             | %           |
| Alkaline Phosphatase     |  0.0019    | Aging       | U/L           | U/L         |
| WBC                      |  0.0554    | Aging       | 10^3/uL       | 10^3/uL     |

**Intercept:** -19.9067
**Chronological age coefficient:** 0.0804

### 6.3 Unit Conversion

The module accepts standard US clinical units and converts internally:

```python
UNIT_CONVERSIONS = {
    "albumin": 10.0,        # g/dL -> g/L (multiply by 10)
    "creatinine": 88.4,     # mg/dL -> umol/L (multiply by 88.4)
    "glucose": 1 / 18.016,  # mg/dL -> mmol/L (divide by 18.016)
}
```

Other biomarkers (lymphocyte %, MCV, RDW, alkaline phosphatase, WBC) use the same units in US and SI systems.

### 6.4 The PhenoAge Formula

**Step 1: Compute the linear predictor (xb)**

```
xb = INTERCEPT + SUM(coefficient_i * SI_value_i) + 0.0804 * chronological_age
```

**Step 2: Compute mortality score via Gompertz model**

```
mortality_score = 1 - exp((MORT_NUMERATOR * exp(xb)) / MORT_DENOMINATOR)
```

Where:
- `MORT_NUMERATOR = -1.51714`  (derived from `-(exp(120 * gamma) - 1)`)
- `MORT_DENOMINATOR = 0.007692696`  (Gompertz shape parameter gamma)

**Step 3: Convert mortality score to biological age**

```
inner = BA_NUMERATOR * ln(1 - mortality_score)
biological_age = (ln(inner) / BA_DENOMINATOR) + BA_INTERCEPT
```

Where:
- `BA_NUMERATOR = -0.0055305`
- `BA_DENOMINATOR = 0.09165`
- `BA_INTERCEPT = 141.50225`

**Step 4: Age acceleration**

```
age_acceleration = biological_age - chronological_age
```

### 6.5 Confidence Intervals

Standard error depends on biomarker completeness:

- All 9 biomarkers available: SE = 4.9 years (from NHANES III validation)
- Fewer than 9 biomarkers: SE = 6.5 years (increased uncertainty)

95% CI: `biological_age +/- 1.96 * SE`

### 6.6 Risk Classification

| Age Acceleration | Risk Level | Meaning                           |
|-----------------|------------|-----------------------------------|
| > +5 years      | HIGH       | Significantly accelerated aging   |
| > +2 years      | MODERATE   | Mildly accelerated aging          |
| -2 to +2 years  | NORMAL     | Aging at expected rate            |
| < -2 years      | LOW        | Aging slower than expected        |

### 6.7 GrimAge Surrogate Estimation

True GrimAge requires DNA methylation data. This module provides a surrogate estimate using plasma proteins that correlate with DNAm GrimAge components (r-squared = 0.72, Hillary et al. 2020, PMID: 32941527).

**Six plasma protein markers:**

| Marker     | Weight | Unit   | Ref Max  |
|------------|--------|--------|----------|
| GDF-15     | 0.15   | pg/mL  | 1,200    |
| Cystatin C | 0.12   | mg/L   | 1.0      |
| PAI-1      | 0.10   | ng/mL  | 43.0     |
| ADM        | 0.11   | pmol/L | 50.0     |
| TIMP-1     | 0.09   | ng/mL  | 250.0    |
| Leptin     | 0.08   | ng/mL  | 15.0     |

**Surrogate formula:**

```
deviation_i = (value_i - ref_max_i) / ref_max_i
weighted_deviation = SUM(weight_i * deviation_i) / SUM(weight_i)
estimated_acceleration = weighted_deviation * 10.0  # empirical scale factor
grimage_score = chronological_age + estimated_acceleration
```

Validation: SE = 5.8 years, from Lothian Birth Cohort 1936 (n=906).

### 6.8 Code Example: Full Calculation

```python
from src.biological_age import BiologicalAgeCalculator

calc = BiologicalAgeCalculator()
result = calc.calculate(
    chronological_age=45,
    biomarkers={
        "albumin": 4.1,             # g/dL
        "creatinine": 0.9,          # mg/dL
        "glucose": 95,              # mg/dL
        "hs_crp": 1.2,             # mg/L (auto-converted to ln_crp)
        "lymphocyte_pct": 30,       # %
        "mcv": 89,                  # fL
        "rdw": 13.5,               # %
        "alkaline_phosphatase": 65, # U/L
        "wbc": 6.5,                # 10^3/uL
        # GrimAge surrogate markers
        "gdf15": 800,              # pg/mL
        "cystatin_c": 0.85,        # mg/L
    },
)
print(f"PhenoAge: {result['biological_age']}")
print(f"Acceleration: {result['age_acceleration']:+.1f} years")
print(f"GrimAge: {result['grimage']['grimage_score']}")
```

---

## Chapter 7: Disease Trajectory Prediction

**File:** `src/disease_trajectory.py` (1,421 lines)

### 7.1 Overview

The `DiseaseTrajectoryAnalyzer` detects pre-symptomatic disease trajectories across 9 disease categories using genotype-stratified biomarker thresholds. It identifies patients on a trajectory toward clinical disease years before conventional diagnosis, enabling early intervention.

### 7.2 The Nine Disease Categories

| Category       | Display Name                    | Key Biomarkers                                      |
|----------------|--------------------------------|-----------------------------------------------------|
| `type2_diabetes`| Type 2 Diabetes               | HbA1c, fasting glucose, fasting insulin, HOMA-IR    |
| `cardiovascular`| Cardiovascular Disease         | Lp(a), LDL-C, ApoB, hs-CRP, TC, HDL-C, TG         |
| `liver`        | Liver Disease (NAFLD/Fibrosis) | ALT, AST, GGT, ferritin, platelets, albumin         |
| `thyroid`      | Thyroid Dysfunction            | TSH, free T4, free T3                               |
| `iron`         | Iron Metabolism Disorder       | Ferritin, transferrin saturation, serum iron, TIBC   |
| `nutritional`  | Nutritional Deficiency         | Omega-3 index, vitamin D, B12, folate, Mg, Zn, Se   |
| `kidney`       | Kidney Disease                 | Creatinine, eGFR, BUN, albumin, cystatin C          |
| `bone_health`  | Bone Health                    | Vitamin D, calcium, PTH, phosphorus                 |
| `cognitive`    | Cognitive Decline              | Homocysteine, B12, folate, hs-CRP, HbA1c            |

### 7.3 Genetic Modifiers

Each disease category includes genetic modifiers that shift risk thresholds:

```python
"type2_diabetes": {
    "genetic_modifiers": {
        "TCF7L2_rs7903146": {"risk_allele": "T", "effect": "beta_cell_dysfunction"},
        "PPARG_rs1801282":  {"risk_allele": "C", "effect": "insulin_sensitivity"},
        "SLC30A8_rs13266634": {"risk_allele": "C", "effect": "zinc_transport"},
        "KCNJ11_rs5219":   {"risk_allele": "T", "effect": "potassium_channel"},
        "GCKR_rs780094":   {"risk_allele": "T", "effect": "glucokinase_regulation"},
    },
}
```

When a patient carries a risk allele, the biomarker thresholds shift -- for example, an HbA1c of 5.7% might be classified as "pre-diabetic" for a TCF7L2 TT carrier but "early metabolic shift" for a CC carrier.

### 7.4 Progression Staging

Each disease has defined stages representing the trajectory from healthy to clinical disease:

```
Type 2 Diabetes: normal -> early_metabolic_shift -> insulin_resistance -> pre_diabetic -> diabetic
Cardiovascular:  optimal -> borderline -> elevated_risk -> high_risk
Liver:           normal -> steatosis_risk -> early_fibrosis -> advanced_fibrosis
Thyroid:         euthyroid -> subclinical -> overt_dysfunction
Iron:            normal -> early_accumulation -> iron_overload
```

### 7.5 Risk Score Formula

The disease trajectory engine computes a composite risk score for each disease category:

1. **Biomarker deviation score**: For each relevant biomarker, calculate deviation from normal range, weighted by clinical importance.
2. **Genetic risk multiplier**: If the patient carries risk alleles, multiply the base risk by a gene-specific factor (typically 1.2x to 2.0x per risk allele).
3. **Age/sex adjustment**: Age and sex modifiers shift thresholds based on epidemiological data.
4. **Composite score**: Weighted combination mapped to risk levels (NORMAL, LOW, MODERATE, HIGH, CRITICAL).

### 7.6 The `analyze_all()` Method

```python
analyzer = DiseaseTrajectoryAnalyzer()
trajectories = analyzer.analyze_all(
    biomarkers={"hba1c": 5.8, "fasting_glucose": 105, "fasting_insulin": 12},
    genotypes={"TCF7L2_rs7903146": "CT"},
    age=45,
    sex="M",
)
# Returns: [
#   {
#     "disease": "type2_diabetes",
#     "risk_level": "MODERATE",
#     "current_stage": "early_metabolic_shift",
#     "current_markers": {"hba1c": 5.8, "fasting_glucose": 105},
#     "genetic_risk_factors": [
#       {"gene": "TCF7L2", "genotype": "CT", "effect": "beta_cell_dysfunction"}
#     ],
#     "years_to_onset_estimate": 8.5,
#     "recommendations": ["Monitor HbA1c every 3 months", "Consider metformin discussion"],
#   },
#   # ... results for other disease categories
# ]
```

### 7.7 Years-to-Onset Estimation

The engine estimates time to clinical onset based on current biomarker levels, rate of change (if longitudinal data available), and genetic risk factors. This is a rough estimate intended to motivate preventive action, not a precise prediction.

---

## Chapter 8: Genotype-Based Reference Ranges

**File:** `src/genotype_adjustment.py` (1,225 lines)

### 8.1 Why Genotype-Adjusted Ranges?

Standard laboratory reference ranges are population averages. Genetic variants can significantly alter what is "normal" for an individual. For example:

- **MTHFR C677T (rs1801133)**: Homozygous TT carriers have 70% reduced enzyme activity, leading to elevated homocysteine. A homocysteine of 12 umol/L is "normal" by standard ranges but may be pathological for a TT carrier.
- **APOE E4**: Carriers have naturally higher LDL-C and respond differently to statin therapy.
- **PNPLA3 I148M**: GG homozygotes have 3x higher risk of NAFLD; their ALT reference range should be tighter.

### 8.2 Core Architecture

The `GenotypeAdjuster` class:

1. Looks up the patient's genotype for each biomarker-gene pair in `GENOTYPE_THRESHOLDS` (from `knowledge.py`).
2. Applies a genotype-specific multiplier or offset to the standard reference range.
3. Returns both the standard and adjusted ranges for comparison.

### 8.3 Ancestry-Specific Adjustments

The `apply_ancestry_adjustments()` method modifies reference ranges based on reported ancestry using data from NHANES III, UK Biobank, and MESA studies. For example:

- eGFR uses the CKD-EPI 2021 equation without race adjustment (PMID: 34554658)
- Vitamin D reference ranges differ by latitude and melanin-mediated synthesis
- Hemoglobin/hematocrit have ancestry-specific normal ranges

### 8.4 Age-Stratified Reference Ranges

Five age brackets with sex-specific ranges:

| Bracket | Age Range | Source Studies                            |
|---------|-----------|------------------------------------------|
| 0-17    | Pediatric | Pediatric guidelines                     |
| 18-39   | Young adult| NHANES III, Framingham                  |
| 40-59   | Middle age | NHANES III, Framingham                  |
| 60-79   | Older adult| KDIGO 2012, ATA 2017, ACC/AHA 2019     |
| 80+     | Elderly    | Geriatric-specific guidelines            |

Example for creatinine:

```python
"creatinine": {
    "18-39": {
        "M": {"low": 0.7, "high": 1.2, "note": "Standard adult male range."},
        "F": {"low": 0.5, "high": 1.0, "note": "Standard adult female range."},
    },
    "60-79": {
        "M": {"low": 0.8, "high": 1.4, "note": "Higher normal; age-related GFR decline."},
        "F": {"low": 0.6, "high": 1.2, "note": "Higher normal; age-related GFR decline."},
    },
}
```

### 8.5 Carrier Screening Integration

For Ashkenazi Jewish patients, the adjuster integrates carrier screening results for compound risk assessment. For example, GBA heterozygous carriers with APOE E4 have a synergistic increase in Parkinson's disease risk.

### 8.6 The `adjust_all()` Method

```python
adjuster = GenotypeAdjuster()
result = adjuster.adjust_all(
    biomarkers={"homocysteine": 12.0, "ldl_c": 145},
    genotypes={"rs1801133": "TT", "APOE": "E3/E4"},
)
# Returns: {
#   "adjustments": [
#     {
#       "biomarker": "homocysteine",
#       "standard_range": {"lower": 5.0, "upper": 15.0},
#       "adjusted_range": {"lower": 5.0, "upper": 10.0},
#       "unit": "umol/L",
#       "gene_display_name": "MTHFR",
#       "genotype_value": "TT",
#       "rationale": "MTHFR 677TT reduces enzyme activity by ~70%; ..."
#     },
#   ]
# }
```

---

## Chapter 9: Clinical Intelligence Modules

Three small, focused modules that provide clinical decision support.

### 9.1 Critical Values Engine

**File:** `src/critical_values.py` (179 lines)

The `CriticalValueEngine` checks biomarker values against life-threatening thresholds that require immediate clinical action. These are distinct from standard reference ranges -- a critical value means "call the physician now."

```python
engine = CriticalValueEngine()
alerts = engine.check({
    "potassium": 6.5,      # Critical high (normal: 3.5-5.0)
    "glucose": 35,         # Critical low (hypoglycemia)
    "sodium": 118,         # Critical low (severe hyponatremia)
})
# Returns: [
#   CriticalValueAlert(biomarker="potassium", value=6.5,
#     threshold_type="HIGH", message="CRITICAL: Potassium 6.5 mEq/L ...")
# ]
```

### 9.2 Discordance Detector

**File:** `src/discordance_detector.py` (299 lines)

The `DiscordanceDetector` identifies contradictions between related biomarkers that suggest a hidden condition or lab error. It implements clinically validated discordance patterns.

**Example discordance patterns:**

| Pattern                    | Biomarkers              | Clinical Implication                    |
|---------------------------|-------------------------|-----------------------------------------|
| LDL/ApoB discordance     | LDL-C low, ApoB high   | Small dense LDL particles; higher risk  |
| Ferritin/iron discordance | Ferritin high, iron low | Inflammation masking iron deficiency    |
| TSH/T4 discordance       | TSH normal, T4 low      | Central hypothyroidism                  |
| AST/ALT ratio            | AST >> ALT              | Alcoholic vs non-alcoholic liver disease|

```python
detector = DiscordanceDetector()
discordances = detector.check({
    "ldl_c": 95,     # Appears normal
    "apob": 130,     # Elevated (discordant with low LDL-C)
    "lpa": 85,       # Elevated Lp(a)
})
# Returns discordance alerts highlighting the LDL/ApoB mismatch
```

### 9.3 Lab Range Interpreter

**File:** `src/lab_range_interpreter.py` (221 lines)

The `LabRangeInterpreter` distinguishes between **standard reference ranges** (what labs report as "normal") and **optimal ranges** (what evidence suggests is ideal for health). Many biomarkers have a significant gap between "not flagged by the lab" and "truly optimal."

**Example:**

| Biomarker  | Standard Range | Optimal Range | Gap                                |
|------------|----------------|---------------|------------------------------------|
| Vitamin D  | 30-100 ng/mL   | 40-60 ng/mL   | 30-39 is "normal" but suboptimal   |
| Ferritin (M)| 12-300 ng/mL  | 40-150 ng/mL  | 12-39 is "normal" but low-optimal  |
| TSH        | 0.45-4.5 mIU/L | 1.0-2.5 mIU/L | 2.5-4.5 is subclinical territory   |

```python
interpreter = LabRangeInterpreter()
discrepancies = interpreter.get_discrepancies(
    biomarkers={"vitamin_d": 32, "tsh": 3.8},
    sex="F",
)
# Returns comparisons showing that both values are within standard range
# but outside optimal range, with interpretation context.
```

---

## Chapter 10: Export System Deep Dive

**Files:** `src/export.py` (1,392 lines) + `src/report_generator.py` (993 lines)

### 10.1 Export Formats

The export system produces five output formats from the same analysis result:

| Format   | Function                        | Use Case                          |
|----------|---------------------------------|-----------------------------------|
| Markdown | `export_markdown()`             | Human-readable reports            |
| JSON     | `export_json()`                 | Machine-readable structured data  |
| PDF      | `export_pdf()`                  | Clinical reports via ReportLab    |
| CSV      | `export_csv()`                  | Spreadsheet analysis              |
| FHIR R4  | `export_fhir_diagnostic_report()` | EHR integration                |

### 10.2 The 12-Section Report

The `ReportGenerator` class produces a structured clinical report:

| Section | Title                            | Content                                    |
|---------|----------------------------------|--------------------------------------------|
| 1       | Biological Age Assessment        | PhenoAge, GrimAge, acceleration, drivers   |
| 2       | Executive Findings               | Top 5 critical/high-priority findings      |
| 3       | Biomarker-Gene Correlation Map   | Which genes affect which biomarkers        |
| 4       | Disease Trajectory Analysis      | Risk for all 9 disease categories          |
| 5       | Pharmacogenomic Profile          | All PGx results with drug recommendations |
| 6       | Nutritional Analysis             | Genotype-aware nutrition assessment        |
| 7       | Interconnected Pathways          | Cross-domain pathway connections           |
| 8       | Prioritized Action Plan          | Ranked interventions by urgency            |
| 9       | Monitoring Schedule              | Follow-up testing timeline                 |
| 10      | Supplement Protocol Summary      | Genotype-guided supplement suggestions     |
| 11      | Clinical Summary for MD          | Concise physician-oriented summary         |
| 12      | References                       | CPIC, PMID citations, data sources         |

### 10.3 PDF Generation via ReportLab

PDF reports use ReportLab's Platypus layout engine:

```python
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter

def export_pdf(query, response_text, evidence=None, analysis=None):
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()

    story = []
    story.append(Paragraph("Biomarker Intelligence Report", styles["Title"]))
    # ... build story elements for each section
    doc.build(story)
    return buffer.getvalue()
```

### 10.4 FHIR R4 DiagnosticReport

The `export_fhir_diagnostic_report()` function produces a FHIR R4 Bundle containing:

- **DiagnosticReport** -- The overall report resource with status, code, and conclusion.
- **Observation** resources -- One per biomarker result, with value, unit, reference range, and interpretation code.
- **Bundle** wrapper -- Transaction bundle for EHR submission.

```python
fhir_bundle = export_fhir_diagnostic_report(
    patient_id="patient-001",
    analysis=analysis_result,
    practitioner_id="dr-smith-001",
)
# Returns: {
#   "resourceType": "Bundle",
#   "type": "transaction",
#   "entry": [
#     {"resource": {"resourceType": "DiagnosticReport", ...}},
#     {"resource": {"resourceType": "Observation", ...}},
#     ...
#   ]
# }
```

### 10.5 Timestamped Filenames

Exported files use UUID-suffixed timestamps to prevent collisions:

```python
generate_filename("pdf")
# -> "biomarker_report_20260311T143025Z_a1b2.pdf"
```

---

## Chapter 11: Testing Strategies

### 11.1 Test Suite Overview

The test suite contains **18 test files** with **709 tests** total. All tests run without external dependencies (Milvus, Claude API) thanks to comprehensive mocking.

**Test distribution by file:**

| File                           | Tests | Focus                                    |
|--------------------------------|------:|------------------------------------------|
| `test_edge_cases.py`           |    69 | Boundary values, malformed inputs, overflow|
| `test_api.py`                  |    59 | FastAPI endpoint testing via TestClient   |
| `test_disease_trajectory.py`   |    48 | All 9 disease categories, staging         |
| `test_export.py`               |    46 | All 5 export formats, content validation  |
| `test_ui.py`                   |    39 | Streamlit component rendering             |
| `test_models.py`               |    39 | Pydantic model validation, serialization  |
| `test_lab_range_interpreter.py`|    37 | Standard vs optimal range comparisons     |
| `test_biological_age.py`       |    30 | PhenoAge formula, GrimAge, edge cases     |
| `test_critical_values.py`      |    28 | Critical threshold alerts                 |
| `test_pharmacogenomics.py`     |    27 | Star allele mapping, drug recommendations |
| `test_genotype_adjustment.py`  |    26 | Genotype and age adjustments              |
| `test_discordance_detector.py` |    25 | Biomarker discordance patterns            |
| `test_collections.py`          |    22 | Schema creation, insert, search           |
| `test_report_generator.py`     |    21 | 12-section report structure               |
| `test_rag_engine.py`           |    21 | RAG pipeline, scoring, prompt building    |
| `test_integration.py`          |    21 | End-to-end agent pipeline                 |
| `test_longitudinal.py`         |    18 | Longitudinal biomarker tracking           |
| `test_agent.py`                |    16 | Agent planning, analysis, synthesis       |

### 11.2 Mock Patterns from conftest.py

The `conftest.py` provides three core fixtures used across all tests:

**Mock Embedder:**

```python
@pytest.fixture
def mock_embedder():
    """Return a mock embedder that produces 384-dim zero vectors."""
    embedder = MagicMock()
    embedder.embed_text.return_value = [0.0] * 384
    return embedder
```

**Mock LLM Client:**

```python
@pytest.fixture
def mock_llm_client():
    """Return a mock LLM client that always responds with 'Mock response'."""
    client = MagicMock()
    client.generate.return_value = "Mock response"
    client.generate_stream.return_value = iter(["Mock ", "response"])
    return client
```

**Mock Collection Manager:**

```python
@pytest.fixture
def mock_collection_manager():
    manager = MagicMock()
    manager.search_all.return_value = {name: [] for name in collection_names}
    manager.get_collection_stats.return_value = {name: 42 for name in collection_names}
    return manager
```

All 14 collections are present in the mock to ensure COLLECTION_CONFIG lookups succeed.

### 11.3 Sample Patient Profile Fixture

```python
@pytest.fixture
def sample_patient():
    return PatientProfile(
        patient_id="TEST-001",
        age=45,
        sex="M",
        biomarkers={
            "albumin": 4.1, "creatinine": 0.9, "glucose": 95,
            "hs_crp": 1.2, "lymphocyte_pct": 30, "mcv": 89,
            "rdw": 13.5, "alkaline_phosphatase": 65, "wbc": 6.5,
        },
        genotypes={"rs1801133": "CT", "APOE": "E3/E4"},
        star_alleles={"CYP2D6": "*1/*4", "CYP2C19": "*1/*2"},
    )
```

### 11.4 Testing Pure Computation Modules

Modules like `biological_age.py`, `pharmacogenomics.py`, `disease_trajectory.py`, and `genotype_adjustment.py` are pure computation -- no I/O, no mocking needed:

```python
def test_phenoage_known_values():
    calc = BiologicalAgeCalculator()
    result = calc.calculate_phenoage(
        chronological_age=50,
        biomarkers={
            "albumin": 4.0, "creatinine": 1.0, "glucose": 100,
            "hs_crp": 2.0, "lymphocyte_pct": 28, "mcv": 90,
            "rdw": 14.0, "alkaline_phosphatase": 70, "wbc": 7.0,
        },
    )
    assert 40 < result["biological_age"] < 70
    assert "mortality_score" in result
    assert len(result["top_aging_drivers"]) <= 5
```

### 11.5 Testing the API

API tests use FastAPI's `TestClient`:

```python
from fastapi.testclient import TestClient
from api.main import app

client = TestClient(app)

def test_health_endpoint():
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert "status" in data
    assert "collections" in data
```

### 11.6 Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=src --cov-report=html

# Run a specific module
pytest tests/test_biological_age.py -v

# Run tests matching a keyword
pytest tests/ -k "phenoage" -v
```

---

## Chapter 12: The Autonomous Agent Pipeline

**File:** `src/agent.py` (610 lines)

### 12.1 Agent Architecture

The `PrecisionBiomarkerAgent` implements the `plan -> analyze -> search -> synthesize -> report` pattern. It wraps the multi-collection RAG engine with four analysis modules and reasoning capabilities.

```
Question + PatientProfile
    |
    v
[Phase 1] analyze_patient()  -- Run all 4 analysis modules
    |   - BiologicalAgeCalculator
    |   - DiseaseTrajectoryAnalyzer
    |   - PharmacogenomicMapper
    |   - GenotypeAdjuster
    |   - CriticalValueEngine
    |   - DiscordanceDetector
    |   - LabRangeInterpreter
    |
    v
[Phase 2] search_plan()  -- Determine search strategy
    |
    v
[Phase 3] rag_engine.retrieve()  -- Multi-collection vector search
    |
    v
[Phase 4] evaluate_evidence()  -- Quality check
    |
    v
[Phase 5] Sub-question expansion (if evidence insufficient)
    |
    v
[Phase 6] _build_enhanced_prompt()  -- Combine evidence + analysis
    |
    v
[Phase 7] LLM generation  -- Claude response
    |
    v
AgentResponse (answer, evidence, analysis, alerts)
```

### 12.2 The SearchPlan Dataclass

```python
@dataclass
class SearchPlan:
    question: str
    identified_topics: List[str] = field(default_factory=list)
    disease_areas: List[str] = field(default_factory=list)
    relevant_modules: List[str] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, domain-specific
    sub_questions: List[str] = field(default_factory=list)
```

### 12.3 Strategy Selection

The agent selects a search strategy based on the question content:

| Strategy         | Condition                                              | Behavior                       |
|------------------|--------------------------------------------------------|-------------------------------|
| `domain-specific`| Single disease area, 0-1 analysis modules              | Focused collection subset     |
| `targeted`       | Specific analysis modules identified                   | Module-guided search          |
| `broad`          | No specific domain or module detected                  | Search all 14 collections     |

### 12.4 Sub-Question Decomposition

Complex questions are decomposed into sub-questions:

- **"Why is X elevated?"** generates:
  - "What genetic variants cause elevated biomarker levels?"
  - "What lifestyle factors contribute to elevated biomarker levels?"
  - "What medications affect biomarker levels?"

- **"Compare X vs Y"** generates:
  - "What are the differences in clinical interpretation?"
  - "What are the genotype-specific considerations?"

- **"What supplements/treatments for X?"** generates:
  - "What are the evidence-based interventions for this condition?"
  - "What genetic factors affect treatment response?"

### 12.5 Evidence Quality Evaluation

```python
def evaluate_evidence(self, evidence: CrossCollectionResult) -> str:
    if evidence.hit_count == 0:
        return "insufficient"
    collections_with_hits = len(evidence.hits_by_collection())
    if collections_with_hits >= 3 and evidence.hit_count >= 10:
        return "sufficient"
    elif collections_with_hits >= 2 and evidence.hit_count >= 5:
        return "partial"
    else:
        return "insufficient"
```

When evidence is "insufficient" and sub-questions exist, the agent runs up to 2 additional retrieval passes with decomposed sub-questions and merges the results.

### 12.6 Critical Alert Extraction

The agent extracts critical alerts from analysis results:

- **Biological age acceleration > 5 years**: "CRITICAL: Biological age acceleration of X years..."
- **Disease trajectory at HIGH/CRITICAL**: "HIGH RISK: cardiovascular trajectory at high level..."
- **DPYD poor/intermediate metabolizer**: "CRITICAL PGx: DPYD -- fluoropyrimidine toxicity risk..."
- **CYP2D6 ultra-rapid**: "PGx ALERT: CYP2D6 -- avoid codeine/tramadol..."
- **CYP2C19 poor/intermediate**: "PGx ALERT: CYP2C19 -- clopidogrel may be ineffective..."
- **Critical value thresholds**: From `CriticalValueEngine`
- **Biomarker discordances**: From `DiscordanceDetector`
- **Optimization opportunities**: From `LabRangeInterpreter`
- **Age-adjusted flags**: From `GenotypeAdjuster.apply_age_adjustments()`

### 12.7 Full Usage Example

```python
from src.agent import PrecisionBiomarkerAgent
from src.models import PatientProfile

agent = PrecisionBiomarkerAgent(rag_engine=engine)

profile = PatientProfile(
    patient_id="PAT-001",
    age=52,
    sex="M",
    biomarkers={
        "albumin": 3.8, "creatinine": 1.1, "glucose": 112,
        "hs_crp": 3.5, "hba1c": 5.9, "ldl_c": 155,
        "apob": 135, "lpa": 85,
    },
    genotypes={
        "TCF7L2_rs7903146": "CT",
        "APOE": "E3/E4",
        "rs1801133": "CT",
    },
    star_alleles={
        "CYP2D6": "*1/*4",
        "CYP2C19": "*1/*2",
    },
)

response = agent.run(
    question="Assess my cardiovascular and metabolic risk profile",
    patient_profile=profile,
)

print(response.answer)
print(f"Critical alerts: {len(response.critical_alerts)}")
print(f"PGx results: {len(response.pgx_results)}")
print(f"Bio age: {response.biological_age.biological_age:.1f}")
```

---

## Chapter 13: Production Deployment

### 13.1 Docker Multi-Stage Build

The Dockerfile uses a two-stage build to minimize image size:

**Stage 1 (builder):** Installs build tools (`gcc`, `g++`) and compiles Python dependencies into a virtual environment at `/opt/venv`.

**Stage 2 (runtime):** Copies only the compiled venv and application source. Runs as non-root user `biomarkeruser`.

```dockerfile
# Stage 1: Build dependencies
FROM python:3.10-slim AS builder
WORKDIR /build
RUN apt-get update && apt-get install -y build-essential gcc g++ ...
COPY requirements.txt .
RUN python -m venv /opt/venv && pip install -r requirements.txt

# Stage 2: Runtime
FROM python:3.10-slim
COPY --from=builder /opt/venv /opt/venv
COPY src/ api/ app/ config/ scripts/ data/ /app/
RUN useradd -r -s /bin/false biomarkeruser
USER biomarkeruser
EXPOSE 8528 8529
```

### 13.2 Compose Topology

The agent runs alongside the HCLS AI Factory services in `docker-compose.dgx-spark.yml`:

```yaml
biomarker-agent:
  build: ./ai_agent_adds/precision_biomarker_agent
  ports:
    - "8528:8528"  # Streamlit UI
    - "8529:8529"  # FastAPI API
  environment:
    - BIOMARKER_MILVUS_HOST=milvus-standalone
    - BIOMARKER_MILVUS_PORT=19530
    - ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
  depends_on:
    - milvus-standalone
    - etcd
    - minio
  healthcheck:
    test: ["CMD", "python", "-c",
           "import urllib.request; urllib.request.urlopen('http://localhost:8528/health')"]
    interval: 30s
    timeout: 10s
    start_period: 60s
    retries: 3
```

### 13.3 Health Checks

The API provides health checks at two levels:

**`GET /health`** -- Returns collection count, total vector count, and agent readiness:

```json
{
  "status": "healthy",
  "collections": 14,
  "total_vectors": 2847,
  "agent_ready": true
}
```

**Docker HEALTHCHECK** -- Uses Python's `urllib` (no `curl` dependency) to probe the Streamlit health endpoint every 30 seconds.

### 13.4 Prometheus Monitoring

The `GET /metrics` endpoint exposes Prometheus-compatible counters:

```
biomarker_api_requests_total 1234
biomarker_api_query_requests_total 567
biomarker_api_analyze_requests_total 89
biomarker_api_errors_total 3
biomarker_collection_vectors{collection="biomarker_reference"} 150
biomarker_collection_vectors{collection="biomarker_genetic_variants"} 320
```

### 13.5 Security Considerations

- **API Key Authentication**: When `BIOMARKER_API_KEY` is set, all endpoints (except `/health` and `/metrics`) require `X-API-Key` header.
- **Request Size Limiting**: Middleware rejects requests exceeding `BIOMARKER_MAX_REQUEST_SIZE_MB` (default 10 MB).
- **CORS**: Restricted to configured origins (default: localhost ports 8080, 8528, 8529).
- **Non-root container**: Runtime user is `biomarkeruser` with no shell access.
- **Input sanitization**: Milvus filter expressions are validated with a safe-character regex (`^[A-Za-z0-9 _\-]+$`) to prevent injection.

### 13.6 Startup Sequence

```
1. Connect to Milvus (host:port from settings)
2. Load SentenceTransformer model (BAAI/bge-small-en-v1.5)
3. Initialize Anthropic Claude client
4. Load knowledge module (static knowledge graph)
5. Initialize analysis modules (BiologicalAgeCalculator, etc.)
6. Build BiomarkerRAGEngine
7. Build PrecisionBiomarkerAgent
8. Store references on app.state for route access
9. Start accepting requests
```

### 13.7 Graceful Shutdown

On SIGTERM/SIGINT, the lifespan context manager disconnects from Milvus:

```python
@asynccontextmanager
async def lifespan(app: FastAPI):
    # ... startup code ...
    yield
    # Shutdown
    if _manager:
        _manager.disconnect()
```

---

## Chapter 14: Future Architecture

### 14.1 Multi-Agent Coordination

The cross-modal event system (`api/routes/events.py`) is the foundation for multi-agent communication:

- **Biomarker -> Imaging Agent**: Elevated Lp(a) triggers coronary calcium scoring recommendation.
- **Biomarker -> CAR-T/Oncology Agent**: DPYD poor metabolizer PGx alert forwarded to oncology pipeline.
- **Imaging Agent -> Biomarker Agent**: Imaging findings trigger biomarker panel recommendations.
- **Biomarker -> Genomics Pipeline**: Unexpected biomarker patterns trigger VCF re-analysis.

Current implementation uses in-memory event stores. Production would use a message bus (NATS, Kafka, or Redis Streams).

### 14.2 Streaming Biomarker Ingestion

Real-time biomarker ingestion from wearables and continuous monitors:

- CGM (continuous glucose monitoring) data -> real-time trajectory updates
- Wearable HRV and resting heart rate -> cardiovascular risk refinement
- Event-driven re-analysis when new data arrives

### 14.3 Fine-Tuned Embeddings

The current BGE-small-en-v1.5 model is general-purpose. Domain-specific fine-tuning opportunities:

- Fine-tune on ClinVar/PharmGKB/CPIC corpus for better biomedical retrieval
- Contrastive learning on biomarker-gene-drug triplets
- Matryoshka Representation Learning for variable-dimension embeddings (128/256/384)

### 14.4 Longitudinal Analysis

Extending the agent to track biomarker trajectories over time:

- Trend detection (improving/worsening/stable) across multiple lab draws
- Velocity-based risk prediction (rate of change matters more than absolute value)
- Intervention effectiveness monitoring (did the supplement protocol work?)

### 14.5 Federated Learning

Privacy-preserving model improvement across institutions:

- Differential privacy for PhenoAge coefficient refinement
- Federated fine-tuning of the embedding model
- Secure aggregation of trajectory risk models

---

## Appendix A: Complete API Reference

### Root Endpoints

#### `GET /`

Returns service info. No authentication required.

**Response:**
```json
{"service": "Biomarker Intelligence Agent", "docs": "/docs", "health": "/health"}
```

#### `GET /health`

**Response (200):**
```json
{
  "status": "healthy",
  "collections": 14,
  "total_vectors": 2847,
  "agent_ready": true
}
```

**Response (503):** Milvus unavailable.

#### `GET /collections`

**Response (200):**
```json
{
  "collections": [
    {"name": "biomarker_reference", "record_count": 150},
    {"name": "biomarker_genetic_variants", "record_count": 320}
  ],
  "total": 14
}
```

#### `GET /knowledge/stats`

**Response (200):**
```json
{
  "disease_domains": 6,
  "total_biomarkers": 45,
  "total_genetic_modifiers": 28,
  "pharmacogenes": 14,
  "pgx_drug_interactions": 35,
  "phenoage_markers": 9,
  "cross_modal_links": 12
}
```

#### `GET /metrics`

Returns Prometheus-formatted plain text with counters and gauges.

### Analysis Endpoints (`/v1`)

#### `POST /v1/analyze`

Full patient analysis (all modules).

**Request:**
```json
{
  "patient_id": "PAT-001",
  "age": 45,
  "sex": "M",
  "biomarkers": {"albumin": 4.1, "creatinine": 0.9, "glucose": 95},
  "genotypes": {"rs1801133": "CT"},
  "star_alleles": {"CYP2D6": "*1/*4"}
}
```

**Response (200):**
```json
{
  "biological_age": {"chronological_age": 45, "biological_age": 43.2, "age_acceleration": -1.8},
  "disease_trajectories": [{"disease": "diabetes", "risk_level": "low", "current_stage": "normal"}],
  "pgx_results": [{"gene": "CYP2D6", "phenotype": "intermediate", "drugs_affected": [...]}],
  "genotype_adjustments": [{"biomarker": "homocysteine", "standard_range": "5-15", "adjusted_range": "5-12"}],
  "critical_alerts": []
}
```

#### `POST /v1/biological-age`

Biological age calculation only.

**Request:**
```json
{
  "age": 50,
  "biomarkers": {
    "albumin": 4.0, "creatinine": 1.0, "glucose": 100,
    "hs_crp": 2.0, "lymphocyte_pct": 28, "mcv": 90,
    "rdw": 14.0, "alkaline_phosphatase": 70, "wbc": 7.0
  }
}
```

**Response (200):**
```json
{
  "chronological_age": 50,
  "biological_age": 48.3,
  "age_acceleration": -1.7,
  "mortality_score": 0.023456,
  "mortality_risk": "NORMAL",
  "confidence_interval": {"lower": 38.7, "upper": 57.9},
  "top_aging_drivers": [...]
}
```

#### `POST /v1/disease-risk`

Disease trajectory analysis.

**Request:**
```json
{
  "age": 45,
  "sex": "M",
  "biomarkers": {"hba1c": 5.8, "fasting_glucose": 105},
  "genotypes": {"TCF7L2_rs7903146": "CT"}
}
```

**Response (200):** List of disease trajectory results across all 9 categories.

#### `POST /v1/pgx`

Pharmacogenomic mapping.

**Request:**
```json
{
  "star_alleles": {"CYP2D6": "*4/*4", "CYP2C19": "*1/*2"},
  "genotypes": {"rs1801133": "CT"}
}
```

**Response (200):**
```json
{
  "gene_results": [
    {
      "gene": "CYP2D6",
      "star_alleles": "*4/*4",
      "phenotype": "Poor Metabolizer",
      "affected_drugs": [
        {"drug": "codeine", "recommendation": "AVOID codeine...", "action": "AVOID", "alert_level": "CRITICAL"}
      ]
    }
  ]
}
```

#### `POST /v1/query`

RAG Q&A query with optional patient profile.

**Request:**
```json
{
  "question": "What does my HbA1c of 5.8% mean?",
  "patient_profile": {
    "patient_id": "PAT-001",
    "age": 45,
    "sex": "M",
    "biomarkers": {"hba1c": 5.8},
    "genotypes": {"TCF7L2_rs7903146": "CT"},
    "star_alleles": {}
  }
}
```

**Response (200):**
```json
{
  "answer": "Based on the evidence...",
  "evidence": {"query": "...", "hits": [...], "total_collections_searched": 14},
  "search_time_ms": 234.5
}
```

#### `GET /v1/health`

V1-specific health check.

### Report Endpoints (`/v1/report`)

#### `POST /v1/report/generate`

Generate a full 12-section patient report.

**Request:** Same as `/v1/analyze`.

**Response (200):**
```json
{
  "report_id": "rpt-a1b2c3d4",
  "generated_at": "2026-03-11T14:30:25Z",
  "markdown": "# Biomarker Intelligence Report\n\n...",
  "analysis_summary": {...}
}
```

#### `GET /v1/report/{report_id}/pdf`

Download a previously generated report as PDF.

**Response (200):** `application/pdf` binary stream.

#### `POST /v1/report/fhir`

Export analysis as FHIR R4 DiagnosticReport Bundle.

**Response (200):** FHIR R4 JSON Bundle.

### Event Endpoints (`/v1/events`)

#### `POST /v1/events/inbound`

Receive cross-modal event from another agent.

**Request:**
```json
{
  "source_agent": "imaging_intelligence_agent",
  "event_type": "imaging_finding",
  "payload": {"finding": "coronary calcification", "severity": "moderate"},
  "patient_id": "PAT-001"
}
```

#### `GET /v1/events/outbound`

Retrieve pending outbound alerts for other agents.

#### `POST /v1/events/alert`

Send a biomarker alert to the platform event bus.

---

## Appendix B: Configuration Reference

All settings use the `BIOMARKER_` prefix and are defined in `config/settings.py` via Pydantic BaseSettings. They can be set via environment variables or `.env` file.

### Path Settings

| Env Var                    | Type | Default                              | Description                    |
|----------------------------|------|--------------------------------------|--------------------------------|
| `BIOMARKER_DATA_DIR`       | Path | `<project_root>/data`                | Data directory                 |
| `BIOMARKER_CACHE_DIR`      | Path | `<project_root>/data/cache`          | Cache directory                |
| `BIOMARKER_REFERENCE_DIR`  | Path | `<project_root>/data/reference`      | Reference data directory       |

### Milvus Settings

| Env Var                    | Type | Default     | Description                   |
|----------------------------|------|-------------|-------------------------------|
| `BIOMARKER_MILVUS_HOST`    | str  | `localhost` | Milvus server hostname        |
| `BIOMARKER_MILVUS_PORT`    | int  | `19530`     | Milvus server port            |
| `BIOMARKER_MILVUS_TIMEOUT_SECONDS` | int | `10` | Milvus operation timeout     |

### Embedding Settings

| Env Var                          | Type | Default                  | Description              |
|----------------------------------|------|--------------------------|--------------------------|
| `BIOMARKER_EMBEDDING_MODEL`      | str  | `BAAI/bge-small-en-v1.5` | Sentence Transformer model|
| `BIOMARKER_EMBEDDING_DIMENSION`  | int  | `384`                    | Embedding vector size    |
| `BIOMARKER_EMBEDDING_BATCH_SIZE` | int  | `32`                     | Batch size for encoding  |

### LLM Settings

| Env Var                          | Type | Default           | Description               |
|----------------------------------|------|-------------------|---------------------------|
| `BIOMARKER_LLM_PROVIDER`        | str  | `anthropic`       | LLM provider name         |
| `BIOMARKER_LLM_MODEL`           | str  | `claude-sonnet-4-6`| Model ID              |
| `BIOMARKER_ANTHROPIC_API_KEY`   | str  | `None`            | Anthropic API key          |
| `BIOMARKER_LLM_MAX_RETRIES`    | int  | `3`               | Max retry attempts         |

### RAG Search Settings

| Env Var                                | Type  | Default | Description                     |
|----------------------------------------|-------|---------|---------------------------------|
| `BIOMARKER_TOP_K_PER_COLLECTION`       | int   | `5`     | Max results per collection      |
| `BIOMARKER_SCORE_THRESHOLD`            | float | `0.4`   | Minimum similarity score        |
| `BIOMARKER_CITATION_HIGH_THRESHOLD`    | float | `0.75`  | Score threshold for "high" relevance |
| `BIOMARKER_CITATION_MEDIUM_THRESHOLD`  | float | `0.60`  | Score threshold for "medium" relevance|

### Collection Weight Settings

| Env Var                                  | Type  | Default | Collection                |
|------------------------------------------|-------|---------|---------------------------|
| `BIOMARKER_WEIGHT_BIOMARKER_REF`         | float | `0.12`  | biomarker_reference       |
| `BIOMARKER_WEIGHT_GENETIC_VARIANTS`      | float | `0.11`  | biomarker_genetic_variants|
| `BIOMARKER_WEIGHT_PGX_RULES`             | float | `0.10`  | biomarker_pgx_rules       |
| `BIOMARKER_WEIGHT_DISEASE_TRAJECTORIES`  | float | `0.10`  | biomarker_disease_trajectories|
| `BIOMARKER_WEIGHT_CLINICAL_EVIDENCE`     | float | `0.09`  | biomarker_clinical_evidence|
| `BIOMARKER_WEIGHT_GENOMIC_EVIDENCE`      | float | `0.08`  | genomic_evidence          |
| `BIOMARKER_WEIGHT_DRUG_INTERACTIONS`     | float | `0.07`  | biomarker_drug_interactions|
| `BIOMARKER_WEIGHT_AGING_MARKERS`         | float | `0.07`  | biomarker_aging_markers   |
| `BIOMARKER_WEIGHT_NUTRITION`             | float | `0.05`  | biomarker_nutrition       |
| `BIOMARKER_WEIGHT_GENOTYPE_ADJUSTMENTS`  | float | `0.05`  | biomarker_genotype_adjustments|
| `BIOMARKER_WEIGHT_MONITORING`            | float | `0.05`  | biomarker_monitoring      |
| `BIOMARKER_WEIGHT_CRITICAL_VALUES`       | float | `0.04`  | biomarker_critical_values |
| `BIOMARKER_WEIGHT_DISCORDANCE_RULES`     | float | `0.04`  | biomarker_discordance_rules|
| `BIOMARKER_WEIGHT_AJ_CARRIER_SCREENING` | float | `0.03`  | biomarker_aj_carrier_screening|

Weights are validated at startup to sum to ~1.0 (tolerance: +/- 0.05).

### Server Settings

| Env Var                          | Type | Default   | Description                   |
|----------------------------------|------|-----------|-------------------------------|
| `BIOMARKER_API_HOST`             | str  | `0.0.0.0` | API bind address              |
| `BIOMARKER_API_PORT`             | int  | `8529`    | API port                      |
| `BIOMARKER_STREAMLIT_PORT`       | int  | `8528`    | Streamlit UI port             |
| `BIOMARKER_METRICS_ENABLED`      | bool | `true`    | Enable Prometheus metrics     |
| `BIOMARKER_CORS_ORIGINS`         | str  | `http://localhost:8080,...` | Comma-separated CORS origins |
| `BIOMARKER_MAX_REQUEST_SIZE_MB`  | int  | `10`      | Max request body size (MB)    |
| `BIOMARKER_REQUEST_TIMEOUT_SECONDS`| int| `60`      | Request timeout               |

### Authentication Settings

| Env Var                | Type | Default | Description                              |
|------------------------|------|---------|------------------------------------------|
| `BIOMARKER_API_KEY`    | str  | `""`    | API key; empty disables auth             |

### Conversation Settings

| Env Var                              | Type | Default | Description                      |
|--------------------------------------|------|---------|----------------------------------|
| `BIOMARKER_MAX_CONVERSATION_CONTEXT` | int  | `3`     | Max conversation turns in memory |

---

## Appendix C: Collection Schema Reference

All collections use `IVF_FLAT` index with `COSINE` metric and `384`-dimensional `FLOAT_VECTOR` embeddings from `BAAI/bge-small-en-v1.5`.

### 1. biomarker_reference

Reference biomarker definitions, ranges, and clinical significance.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique biomarker identifier              |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `name`               | VARCHAR       | 100       | Biomarker display name                   |
| `unit`               | VARCHAR       | 20        | Measurement unit (e.g., mg/dL)           |
| `category`           | VARCHAR       | 30        | CBC, CMP, Lipids, Thyroid, etc.          |
| `ref_range_min`      | FLOAT         | --        | Standard reference range lower bound     |
| `ref_range_max`      | FLOAT         | --        | Standard reference range upper bound     |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |
| `clinical_significance`| VARCHAR     | 2000      | Clinical interpretation                  |
| `epigenetic_clock`   | VARCHAR       | 50        | PhenoAge/GrimAge coefficient if applicable|
| `genetic_modifiers`  | VARCHAR       | 500       | Comma-separated modifier genes           |

### 2. biomarker_genetic_variants

Genetic variants affecting biomarker levels and disease risk.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique variant identifier                |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `gene`               | VARCHAR       | 50        | Gene symbol (e.g., MTHFR)               |
| `rs_id`              | VARCHAR       | 20        | dbSNP rsID (e.g., rs1801133)            |
| `risk_allele`        | VARCHAR       | 20        | Risk allele                              |
| `protective_allele`  | VARCHAR       | 5         | Protective allele                        |
| `effect_size`        | VARCHAR       | 250       | Effect size description                  |
| `mechanism`          | VARCHAR       | 2000      | Molecular mechanism                      |
| `disease_associations`| VARCHAR      | 1000      | Comma-separated disease associations     |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 3. biomarker_pgx_rules

Pharmacogenomic dosing rules following CPIC guidelines.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique PGx rule identifier               |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `gene`               | VARCHAR       | 50        | Pharmacogene (e.g., CYP2D6)             |
| `star_alleles`       | VARCHAR       | 100       | Star allele combination (e.g., *1/*2)    |
| `drug`               | VARCHAR       | 100       | Drug name                                |
| `phenotype`          | VARCHAR       | 30        | Metabolizer phenotype                    |
| `cpic_level`         | VARCHAR       | 10        | CPIC evidence level (1A, 1B, 2A, etc.)  |
| `recommendation`     | VARCHAR       | 2000      | Clinical recommendation text             |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 4. biomarker_disease_trajectories

Disease progression trajectory definitions and staging criteria.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique trajectory identifier             |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `disease`            | VARCHAR       | 50        | Disease category                         |
| `disease_area`       | VARCHAR       | 50        | Disease area for filtering               |
| `stage`              | VARCHAR       | 50        | Progression stage                        |
| `biomarker_pattern`  | VARCHAR       | 2000      | Biomarker criteria for this stage        |
| `genetic_modifiers`  | VARCHAR       | 500       | Genetic modifiers affecting trajectory   |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 5. biomarker_clinical_evidence

Published clinical evidence with PubMed linkage.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique evidence identifier               |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `title`              | VARCHAR       | 500       | Publication title                        |
| `authors`            | VARCHAR       | 500       | Author list                              |
| `year`               | INT64         | --        | Publication year (used for date filters) |
| `pmid`               | VARCHAR       | 20        | PubMed ID                                |
| `disease_area`       | VARCHAR       | 50        | Disease area for filtering               |
| `evidence_level`     | VARCHAR       | 20        | Evidence level classification            |
| `text_chunk`         | VARCHAR       | 3000      | Abstract/summary for embedding           |
| `text_summary`       | VARCHAR       | 2000      | Concise summary                          |

### 6. biomarker_nutrition

Genotype-aware nutritional guidance.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique guideline identifier              |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `nutrient`           | VARCHAR       | 100       | Nutrient name                            |
| `gene`               | VARCHAR       | 50        | Relevant gene                            |
| `genotype`           | VARCHAR       | 20        | Genotype that modifies recommendation    |
| `recommendation`     | VARCHAR       | 2000      | Nutritional recommendation               |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 7. biomarker_drug_interactions

Gene-drug interaction records.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique interaction identifier            |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `drug_name`          | VARCHAR       | 100       | Drug name                                |
| `gene`               | VARCHAR       | 50        | Interacting gene                         |
| `interaction_type`   | VARCHAR       | 50        | Type of interaction                      |
| `severity`           | VARCHAR       | 20        | Severity level                           |
| `recommendation`     | VARCHAR       | 2000      | Clinical recommendation                  |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 8. biomarker_aging_markers

Epigenetic aging clock markers and correlations.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique marker identifier                 |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `marker_name`        | VARCHAR       | 100       | Aging marker name                        |
| `clock_type`         | VARCHAR       | 50        | PhenoAge, GrimAge, etc.                  |
| `coefficient`        | FLOAT         | --        | Clock coefficient value                  |
| `direction`          | VARCHAR       | 20        | Aging or protective                      |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 9. biomarker_genotype_adjustments

Genotype-based reference range adjustments.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique adjustment identifier             |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `biomarker`          | VARCHAR       | 100       | Biomarker being adjusted                 |
| `gene`               | VARCHAR       | 50        | Gene causing adjustment                  |
| `genotype`           | VARCHAR       | 20        | Specific genotype                        |
| `standard_range`     | VARCHAR       | 50        | Standard reference range                 |
| `adjusted_range`     | VARCHAR       | 50        | Genotype-adjusted range                  |
| `rationale`          | VARCHAR       | 2000      | Clinical rationale                       |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 10. biomarker_monitoring

Condition-specific monitoring protocols.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique protocol identifier               |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `condition`          | VARCHAR       | 100       | Condition being monitored                |
| `biomarker`          | VARCHAR       | 100       | Biomarker to monitor                     |
| `frequency`          | VARCHAR       | 50        | Monitoring frequency                     |
| `rationale`          | VARCHAR       | 2000      | Why this monitoring is needed            |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 11. biomarker_critical_values

Critical value thresholds requiring immediate clinical action.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique critical value identifier         |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `biomarker`          | VARCHAR       | 100       | Biomarker name                           |
| `threshold_high`     | FLOAT         | --        | Critical high threshold                  |
| `threshold_low`      | FLOAT         | --        | Critical low threshold                   |
| `unit`               | VARCHAR       | 20        | Unit of measurement                      |
| `action`             | VARCHAR       | 2000      | Required clinical action                 |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 12. biomarker_discordance_rules

Cross-biomarker discordance detection rules.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique rule identifier                   |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `biomarker_a`        | VARCHAR       | 100       | First biomarker in pair                  |
| `biomarker_b`        | VARCHAR       | 100       | Second biomarker in pair                 |
| `pattern`            | VARCHAR       | 500       | Expected vs discordant pattern           |
| `clinical_meaning`   | VARCHAR       | 2000      | Clinical interpretation                  |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 13. biomarker_aj_carrier_screening

Ashkenazi Jewish genetic carrier screening panel.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Unique screening entry identifier        |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `gene`               | VARCHAR       | 50        | Gene name (BRCA1, HEXA, GBA, etc.)      |
| `condition`          | VARCHAR       | 200       | Associated condition                     |
| `carrier_frequency`  | VARCHAR       | 50        | Population carrier frequency             |
| `inheritance`        | VARCHAR       | 50        | Inheritance pattern                      |
| `compound_risks`     | VARCHAR       | 1000      | Compound risk interactions               |
| `text_chunk`         | VARCHAR       | 3000      | Text chunk used for embedding            |

### 14. genomic_evidence (read-only, shared)

Shared genomic variant evidence collection from the VCF-derived pipeline. Read-only for the biomarker agent; written by the genomics pipeline.

| Field                | Type          | Max Length | Description                              |
|----------------------|---------------|-----------|------------------------------------------|
| `id`                 | VARCHAR (PK)  | 100       | Variant identifier                       |
| `embedding`          | FLOAT_VECTOR  | dim=384   | BGE-small-en-v1.5 text embedding         |
| `chrom`              | VARCHAR       | 10        | Chromosome                               |
| `pos`                | INT64         | --        | Genomic position                         |
| `ref`                | VARCHAR       | 500       | Reference allele                         |
| `alt`                | VARCHAR       | 500       | Alternate allele                         |
| `gene`               | VARCHAR       | 50        | Gene symbol                              |
| `consequence`        | VARCHAR       | 100       | Variant consequence (missense, etc.)     |
| `clinvar_significance`| VARCHAR      | 100       | ClinVar clinical significance            |
| `text_chunk`         | VARCHAR       | 3000      | Text summary for embedding               |

---

*End of Learning Guide (Advanced) -- Precision Biomarker Intelligence Agent*
*Total codebase: 12,628 lines source + 8,772 lines tests = 21,400 lines across 36 files.*
