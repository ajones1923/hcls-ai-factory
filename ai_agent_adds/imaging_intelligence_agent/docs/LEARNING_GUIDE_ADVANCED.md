# Imaging Intelligence Agent -- Advanced Learning Guide

**Version:** 1.0
**Date:** March 2026
**Author:** Adam Jones
**Audience:** Senior engineers, ML researchers, and platform architects extending the Imaging Intelligence Agent

---

## Prerequisites

### What You Should Know

Before working through this guide you should be comfortable with:

- **Python 3.10+** -- async/await, dataclasses, type hints, ABC pattern
- **Pydantic v2** -- BaseModel, Field validators, model_dump/model_dump_json
- **Vector databases** -- embedding similarity search, IVF indexes, distance metrics
- **Medical imaging basics** -- DICOM, CT windowing, MRI sequences, modality codes
- **FastAPI** -- request/response models, routers, middleware, dependency injection
- **Docker Compose** -- multi-service orchestration, volumes, networks, health checks
- **pytest** -- fixtures, mocking, parametrize, conftest patterns

### Codebase Map

```
agent/
  api/
    main.py                         # FastAPI server (port 8524)
    routes/
      events.py                     # DICOM webhook + workflow routing
      meta_agent.py                 # Multi-step agent endpoint
      nim.py                        # NIM service status/invoke
      preview.py                    # Preview generation endpoints
      protocol.py                   # Protocol optimization
      dose.py                       # Dose comparison endpoints
      reports.py                    # Export endpoints (MD/JSON/PDF/FHIR)
      workflows.py                  # Workflow execution endpoints
      demo_cases.py                 # Pre-built demo case runner
  app/                              # Streamlit chat UI (port 8525)
  config/
    settings.py                     # ImagingSettings (50+ env vars)
    ohif-config.js                  # OHIF Viewer configuration
  src/
    agent.py                        # ImagingIntelligenceAgent orchestrator
    collections.py                  # ImagingCollectionManager (Milvus)
    cross_modal.py                  # Imaging->Genomics trigger (401 lines)
    demo_cases.py                   # Demo case definitions
    export.py                       # 4-format export (691 lines)
    knowledge.py                    # Domain knowledge graph (1843 lines)
    models.py                       # 37 enums + 20 Pydantic models (728 lines)
    protocol_optimizer.py           # AI protocol optimization
    query_expansion.py              # 12 expansion maps (186 lines)
    rag_engine.py                   # ImagingRAGEngine (692 lines)
    scheduler.py                    # APScheduler ingest scheduler
    imaging/
      preview_generator.py          # MP4/GIF slice animations
    ingest/
      base.py                       # BaseIngestPipeline ABC
      literature_parser.py          # PubMed ingest
      clinical_trials_parser.py     # ClinicalTrials.gov ingest
      finding_parser.py             # Imaging finding templates
      protocol_parser.py            # Acquisition protocol ingest
      device_parser.py              # FDA-cleared device ingest
      anatomy_parser.py             # Anatomical structure ingest
      benchmark_parser.py           # Model benchmark ingest
      guideline_parser.py           # Clinical guideline ingest
      report_template_parser.py     # Structured report templates
      dataset_parser.py             # Public dataset registry
      dicom_watcher.py              # Orthanc change-feed watcher
    nim/
      base.py                       # BaseNIMClient (83 lines)
      llm_client.py                 # LlamaLLMClient (439 lines)
      vista3d_client.py             # VISTA3DClient (215 lines)
      maisi_client.py               # MAISIClient (136 lines)
      vilam3_client.py              # VILAM3Client (537 lines)
      service_manager.py            # NIMServiceManager
    workflows/
      base.py                       # BaseImagingWorkflow (83 lines)
      ct_head_hemorrhage.py         # ICH detection (488 lines)
      ct_chest_lung_nodule.py       # Lung nodule tracking (598 lines)
      ct_coronary_angiography.py    # CAD-RADS scoring (357 lines)
      cxr_rapid_findings.py         # CXR classification (658 lines)
      mri_brain_ms_lesion.py        # MS lesion analysis (645 lines)
      mri_prostate_pirads.py        # PI-RADS scoring (280 lines)
  tests/
    conftest.py                     # 10+ fixtures
    test_models.py                  # Data model validation
    test_export.py                  # Export format tests
    test_rag_engine.py              # RAG engine tests
    test_workflows.py               # Workflow tests
    test_agent.py                   # Agent orchestrator tests
    test_knowledge.py               # Knowledge graph tests
    test_nim_clients.py             # NIM client tests
    test_dicom_ingestion.py         # DICOM watcher tests
    test_cross_modal.py             # Cross-modal trigger tests
    test_query_expansion.py         # Query expansion tests
    test_preview_generator.py       # Preview generation tests
    test_preview_api.py             # Preview API endpoint tests
  flare/                            # NVIDIA FLARE federated learning
    job_configs/
      cxr_classification/           # Chest X-ray classification FL job
      ct_segmentation/              # CT organ segmentation FL job
      lung_nodule_detection/        # Lung nodule detection FL job
  docker-compose.yml                # Full stack (11 services)
  docker-compose.lite.yml           # Lite stack (6 services, no GPU)
  Dockerfile                        # Multi-stage build
```

---

## Chapter 1: Deep Dive into the RAG Engine

The `ImagingRAGEngine` class in `src/rag_engine.py` (692 lines) is the central retrieval-augmented generation component. It coordinates multi-collection vector search, knowledge graph augmentation, query expansion, comparative analysis, and LLM synthesis.

### 1.1 Class Architecture

```python
class ImagingRAGEngine:
    def __init__(self, collection_manager, embedder, llm_client, nim_service_manager=None):
        self.collection_manager = collection_manager  # ImagingCollectionManager
        self.embedder = embedder                      # SentenceTransformer (BGE-small)
        self.llm_client = llm_client                  # LlamaLLMClient
        self.nim_manager = nim_service_manager         # Optional NIMServiceManager
        self.system_prompt = SYSTEM_PROMPT
```

The engine accepts four dependencies via constructor injection, making every component mockable for testing.

### 1.2 The retrieve() Method

`retrieve()` is the core search pipeline. Its execution flow:

1. **Query expansion** -- calls `expand_query(query)` to add synonyms from 12 domain maps
2. **Embedding** -- encodes the expanded query via `_embed_query()` (BGE-small-en-v1.5, 384-dim, L2-normalized)
3. **Knowledge context** -- checks query against `IMAGING_PATHOLOGIES`, `IMAGING_MODALITIES`, and `IMAGING_ANATOMY` dictionaries for supplemental context
4. **Collection iteration** -- searches each collection in `COLLECTION_CONFIG` (or a filtered subset)
5. **Filter construction** -- builds per-collection Milvus filter expressions from modality, body_region, and year parameters
6. **Weighted scoring** -- multiplies each hit's raw similarity score by the collection weight
7. **Sorting** -- sorts all hits by weighted score descending
8. **Result assembly** -- returns a `CrossCollectionResult` with hits, knowledge context, timing, and collection count

```python
def retrieve(self, query, top_k_per_collection=5, collections_filter=None,
             year_min=None, year_max=None, modality_filter=None,
             body_region_filter=None) -> CrossCollectionResult:
```

**Key detail:** The `top_k_per_collection` parameter defaults to 5, meaning a full search across all 11 collections can return up to 55 candidate hits before weighted ranking.

### 1.3 Score Weighting (COLLECTION_CONFIG)

Each collection has a weight that reflects its relevance priority:

| Collection | Label | Weight | Has Modality | Year Field |
|---|---|---|---|---|
| `imaging_literature` | Literature | 0.18 | Yes | `year` |
| `imaging_findings` | Finding | 0.15 | Yes | -- |
| `imaging_trials` | Trial | 0.12 | Yes | `start_year` |
| `imaging_guidelines` | Guideline | 0.10 | Yes | `year` |
| `imaging_protocols` | Protocol | 0.08 | Yes | -- |
| `imaging_devices` | Device | 0.08 | Yes | -- |
| `imaging_benchmarks` | Benchmark | 0.08 | Yes | -- |
| `imaging_anatomy` | Anatomy | 0.06 | No | -- |
| `imaging_datasets` | Dataset | 0.06 | Yes | -- |
| `imaging_report_templates` | ReportTemplate | 0.05 | Yes | -- |
| `genomic_evidence` | Genomic | 0.04 | No | -- |

**Weights sum to 1.00.** Literature and Findings dominate because they contain the most clinically actionable evidence. The genomic collection has the lowest weight (0.04) because it provides cross-modal enrichment context rather than primary imaging evidence.

The weighted score formula is:

```
final_score = raw_cosine_similarity * collection_weight
```

### 1.4 Citation Scoring

Citation quality thresholds are configured in `config/settings.py`:

```python
CITATION_HIGH_THRESHOLD: float = 0.75
CITATION_MEDIUM_THRESHOLD: float = 0.60
```

Hits with `score >= 0.75` are high-confidence citations. Those between 0.60 and 0.75 are medium-confidence. Below 0.60 are low-confidence.

### 1.5 System Prompt

The `SYSTEM_PROMPT` constant (lines 222-241) defines the LLM's persona across 11 knowledge domains:

1. CT Analysis
2. MRI Interpretation
3. Chest X-ray
4. Imaging AI Models (MONAI, VISTA-3D, nnU-Net, SwinUNETR, DenseNet)
5. Clinical Guidelines (ACR, Lung-RADS, BI-RADS, TI-RADS, LI-RADS)
6. Imaging Protocols
7. FDA-Cleared Devices (510(k), De Novo)
8. Radiology Reporting (RadLex, DICOM SR)
9. Public Datasets (RSNA, TCIA, NIH, LIDC-IDRI, BraTS, CheXpert, MIMIC-CXR)
10. Quantitative Imaging (RECIST, volume doubling time)
11. NVIDIA NIMs (VISTA-3D, MAISI, VILA-M3, Llama3)

The prompt enforces two critical rules: (a) always cite evidence from context, and (b) include a research-use-only disclaimer.

### 1.6 Prompt Construction (_build_prompt)

The `_build_prompt()` method assembles the LLM input as a two-message list:

```python
[
    {"role": "system", "content": SYSTEM_PROMPT},
    {"role": "user", "content": "<domain knowledge>\n<evidence>\n<conversation>\n<question>"},
]
```

The user content is built in sections:
- **Domain Knowledge** -- from the knowledge graph (pathology/modality/anatomy)
- **Retrieved Evidence** -- top 20 hits with `[Label] (score: 0.XXX) text` format
- **Conversation Context** -- prior exchanges for multi-turn sessions
- **Question** -- the raw user query

### 1.7 Query Entry Points

Two public methods serve as entry points:

- `query(question, conversation_context, **kwargs) -> str` -- synchronous, full RAG pipeline
- `query_stream(question, conversation_context, **kwargs) -> Generator[str]` -- streaming tokens

Both automatically detect comparative queries via `_is_comparative()` and route to the comparative pipeline when matched.

---

## Chapter 2: Vector Search Internals

### 2.1 Index Type: IVF_FLAT

The Imaging Agent uses `IVF_FLAT` (Inverted File with Flat quantization) defined in `src/collections.py`:

```python
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}
```

IVF_FLAT partitions vectors into 1024 Voronoi cells (clusters). At search time, only a subset of clusters is scanned, controlled by `nprobe`.

**Why IVF_FLAT over HNSW?** For the Imaging Agent's collection sizes (hundreds to low thousands of records per collection), IVF_FLAT offers lower memory overhead than HNSW's graph structure while maintaining exact distance computation within probed clusters. HNSW would be preferred at scale (100K+ vectors per collection).

### 2.2 Distance Metric: COSINE

```python
SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}
```

**COSINE vs L2 vs IP comparison:**

| Metric | Formula | Range | Best For |
|---|---|---|---|
| COSINE | 1 - cos(a,b) | [0, 2] | Normalized embeddings (BGE default) |
| L2 (Euclidean) | \|\|a-b\|\|^2 | [0, inf) | Unnormalized, magnitude matters |
| IP (Inner Product) | a . b | (-inf, inf) | Already-normalized (equivalent to COSINE) |

COSINE is the correct choice because BGE-small-en-v1.5 produces L2-normalized embeddings (via `normalize_embeddings=True`). For normalized vectors, COSINE and IP produce equivalent rankings, but COSINE scores are bounded [0, 1] after Milvus's internal conversion, making threshold comparisons (e.g., `SCORE_THRESHOLD = 0.40`) more interpretable.

### 2.3 nprobe = 16

With `nlist=1024` clusters and `nprobe=16`, each search scans 16/1024 = 1.56% of the index. This provides a good recall/latency tradeoff:

| nprobe | Clusters Scanned | Recall (approx.) | Relative Latency |
|---|---|---|---|
| 1 | 0.1% | ~60% | 1x |
| 8 | 0.8% | ~90% | 4x |
| **16** | **1.56%** | **~95%** | **8x** |
| 64 | 6.25% | ~99% | 32x |
| 1024 | 100% | 100% (exact) | 512x |

For collections under 10K vectors, nprobe=16 effectively achieves near-exact recall because many clusters contain fewer than `top_k` vectors.

### 2.4 BGE Embedding Model

```python
EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
EMBEDDING_DIMENSION: int = 384
EMBEDDING_BATCH_SIZE: int = 32
```

BGE-small-en-v1.5 characteristics:
- **Parameters:** 33M (MiniLM backbone)
- **Dimensions:** 384
- **Max sequence length:** 512 tokens
- **MTEB score:** 63.55 (top-tier for its size class)
- **License:** MIT

The `_embed_query()` method normalizes embeddings:

```python
def _embed_query(self, text: str) -> List[float]:
    return self.embedder.encode(text, normalize_embeddings=True).tolist()
```

### 2.5 Collection Schemas

Each collection defines a `FLOAT_VECTOR` field with `dim=384`. The primary key is always a VARCHAR `id` field. Example from `imaging_literature`:

```python
EMBEDDING_DIM = 384  # BGE-small-en-v1.5

LITERATURE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=50),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="title", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_chunk", dtype=DataType.VARCHAR, max_length=3000),
    # ... additional metadata fields
]
```

---

## Chapter 3: Adding a New Collection

This worked example adds an `imaging_conferences` collection for storing conference proceedings and presentations.

### Step 1: Define the Pydantic Model

Add to `src/models.py`:

```python
class ConferenceRecord(BaseModel):
    """Conference proceeding / presentation abstract."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=3000)
    title: str = Field("", max_length=500)
    conference_name: str = Field("", max_length=200)  # e.g., RSNA, MICCAI
    year: int = Field(0, ge=2000, le=2030)
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    presentation_type: str = Field("", max_length=50)  # poster, oral, keynote
    authors: str = Field("", max_length=500)

    def to_embedding_text(self) -> str:
        parts = [self.title, self.text_summary]
        if self.conference_name:
            parts.append(f"Conference: {self.conference_name}")
        if self.modality:
            parts.append(f"Modality: {self.modality.value}")
        return " ".join(parts)
```

### Step 2: Define the Milvus Schema

Add to `src/collections.py`:

```python
CONFERENCE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="title", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="conference_name", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="year", dtype=DataType.INT64),
    FieldSchema(name="modality", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="body_region", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="presentation_type", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="authors", dtype=DataType.VARCHAR, max_length=500),
]
```

### Step 3: Register the Schema

In the `ImagingCollectionManager.create_all_collections()` method, add:

```python
self._create_collection("imaging_conferences", CONFERENCE_FIELDS)
```

### Step 4: Add to COLLECTION_CONFIG

In `src/rag_engine.py`, add to `COLLECTION_CONFIG`:

```python
"imaging_conferences": {
    "weight": 0.05,
    "label": "Conference",
    "has_modality": True,
    "year_field": "year",
},
```

**Important:** Adjust other weights so they still sum to approximately 1.0.

### Step 5: Add to conftest.py

Update `IMAGING_COLLECTION_NAMES` in `tests/conftest.py`:

```python
IMAGING_COLLECTION_NAMES = [
    # ... existing collections ...
    "imaging_conferences",
]
```

### Step 6: Create the Ingest Parser

Create `src/ingest/conference_parser.py`:

```python
from src.ingest.base import BaseIngestPipeline
from src.models import ConferenceRecord

class ConferenceIngestPipeline(BaseIngestPipeline):
    COLLECTION_NAME = "imaging_conferences"

    def fetch(self, **kwargs):
        # Implement API call or file read
        ...

    def parse(self, raw_data):
        records = []
        for item in raw_data:
            records.append(ConferenceRecord(
                id=item["id"],
                text_summary=item["abstract"],
                title=item["title"],
                conference_name=item["conference"],
                year=item["year"],
            ))
        return records
```

### Step 7: Add Query Expansion Terms

In `src/query_expansion.py`, add a new map or extend an existing one:

```python
CONFERENCE_EXPANSION = {
    "rsna": {"rsna annual meeting", "radiological society", "rsna 2025"},
    "miccai": {"miccai conference", "medical image computing"},
    "isbi": {"isbi conference", "biomedical imaging symposium"},
}
```

Add to `ALL_EXPANSION_MAPS`.

### Step 8: Add Filter Support

In `ImagingRAGEngine.retrieve()`, add `imaging_conferences` to `_has_body_region`:

```python
_has_body_region = {
    # ... existing collections ...
    "imaging_conferences",
}
```

### Step 9: Write Tests

Create `tests/test_conference_ingest.py` with at minimum:
- Model validation tests (valid/invalid records)
- Parse method tests with sample data
- Integration test with mock collection manager

### Step 10: Update Settings

Add collection name to `config/settings.py`:

```python
COLLECTION_CONFERENCES: str = "imaging_conferences"
```

---

## Chapter 4: Building a Custom Ingest Pipeline

### 4.1 The BaseIngestPipeline Pattern

The abstract base class in `src/ingest/base.py` defines a three-step pipeline:

```python
class BaseIngestPipeline(ABC):
    COLLECTION_NAME: str = ""

    def __init__(self, collection_manager, embedder):
        self.collection_manager = collection_manager
        self.embedder = embedder

    @abstractmethod
    def fetch(self, **kwargs) -> Any:
        """Retrieve raw data from external source."""
        ...

    @abstractmethod
    def parse(self, raw_data: Any) -> List[BaseModel]:
        """Convert raw data into validated Pydantic models."""
        ...

    def embed_and_store(self, records, collection_name=None, batch_size=32) -> int:
        """Embed text and insert into Milvus in batches."""
        ...

    def run(self, collection_name=None, batch_size=32, **fetch_kwargs) -> int:
        """Orchestrate fetch -> parse -> embed_and_store."""
        ...
```

### 4.2 The fetch -> parse -> embed_and_store Flow

**fetch()** -- retrieves raw data. Each parser implements this differently:
- `literature_parser.py` calls the PubMed E-utilities API
- `clinical_trials_parser.py` calls the ClinicalTrials.gov v2 API
- `device_parser.py` reads FDA device listing files
- `anatomy_parser.py` reads local JSONL anatomy reference files

**parse()** -- converts raw API/file data into validated Pydantic model instances. Each record must implement `to_embedding_text() -> str`, which produces the string that gets embedded.

**embed_and_store()** -- the base class handles this uniformly:

```python
def embed_and_store(self, records, collection_name=None, batch_size=32):
    for i in range(0, len(records), batch_size):
        batch = records[i : i + batch_size]
        texts = [r.to_embedding_text() for r in batch]
        embeddings = self.embedder.encode(texts, normalize_embeddings=True).tolist()

        insert_data = []
        for rec, emb in zip(batch, embeddings):
            row = rec.model_dump()
            for k, v in row.items():
                if hasattr(v, "value"):
                    row[k] = v.value  # Convert enums to strings
            row["embedding"] = emb
            insert_data.append(row)

        self.collection_manager.insert_batch(coll_name, insert_data, batch_size=batch_size)
```

### 4.3 Batch Embedding

Embedding is done in batches of 32 (configurable via `EMBEDDING_BATCH_SIZE`). This balances GPU memory usage against throughput. For BGE-small-en-v1.5 on a DGX Spark:

| Batch Size | Throughput (records/sec) | GPU Memory |
|---|---|---|
| 1 | ~50 | ~200 MB |
| 16 | ~400 | ~400 MB |
| 32 | ~650 | ~600 MB |
| 64 | ~800 | ~1.0 GB |

### 4.4 Worked Example: FDA Device Parser

```python
class DeviceIngestPipeline(BaseIngestPipeline):
    COLLECTION_NAME = "imaging_devices"

    def fetch(self, data_dir=None, **kwargs):
        # Read FDA device JSON files from reference directory
        path = Path(data_dir or self.default_data_dir)
        devices = []
        for f in path.glob("*.json"):
            with open(f) as fh:
                devices.extend(json.load(fh))
        return devices

    def parse(self, raw_data):
        records = []
        for item in raw_data:
            try:
                records.append(ImagingDevice(
                    id=item["k_number"],
                    text_summary=item["device_description"],
                    device_name=item["device_name"],
                    manufacturer=item["applicant"],
                    regulatory_status=DeviceRegulatory.CLEARED_510K,
                    clearance_date=item.get("decision_date", ""),
                    modality=self._map_modality(item),
                    intended_use=item.get("intended_use", ""),
                ))
            except ValidationError as e:
                logger.warning(f"Skipping invalid device {item.get('k_number')}: {e}")
        return records
```

---

## Chapter 5: Extending the Knowledge Graph

### 5.1 Three Domains

The knowledge graph in `src/knowledge.py` (1843 lines) contains three domain dictionaries:

| Domain | Variable | Entries | Key Fields |
|---|---|---|---|
| Pathologies | `IMAGING_PATHOLOGIES` | ~25 | ICD-10, modalities, CT/MRI characteristics, severity criteria, AI models, guidelines |
| Modalities | `IMAGING_MODALITIES` | ~9 | Physics, protocols, strengths, limitations, AI applications, DICOM codes |
| Anatomy | `IMAGING_ANATOMY` | ~21 | Body region, systems, SNOMED codes, FMA IDs, VISTA-3D labels, common pathologies |

### 5.2 Adding a New Pathology

To add a new pathology (e.g., pulmonary embolism):

```python
IMAGING_PATHOLOGIES["pulmonary_embolism"] = {
    "icd10": "I26.99",
    "display_name": "Pulmonary Embolism",
    "modalities": ["ct"],
    "body_region": "chest",
    "subtypes": ["acute", "chronic", "saddle", "segmental", "subsegmental"],
    "ct_characteristics": (
        "Filling defect in pulmonary artery on CTPA. "
        "Acute: central low-attenuation defect with vessel expansion. "
        "Chronic: eccentric mural thrombus with vessel narrowing."
    ),
    "severity_criteria": {
        "critical": "Saddle PE, RV strain, hemodynamic instability",
        "urgent": "Lobar PE, RV/LV ratio > 1.0",
        "significant": "Segmental PE without RV strain",
        "routine": "Subsegmental PE, no hemodynamic compromise",
    },
    "ai_models": ["MONAI PE Detection", "AIDOC PE Triage"],
    "related_guidelines": ["Fleischner Society CTPA", "AHA PE Guidelines"],
    "genomic_links": ["Factor V Leiden", "Prothrombin G20210A"],
}
```

### 5.3 Adding a New Modality

```python
IMAGING_MODALITIES["nuclear"] = {
    "display_name": "Nuclear Medicine",
    "physics": "Gamma ray detection from administered radiotracers",
    "dicom_code": "NM",
    "typical_studies": ["bone scan", "thyroid scan", "renal scan", "MIBG"],
    # ...
}
```

### 5.4 Adding a New Anatomy Entry

```python
IMAGING_ANATOMY["aorta"] = {
    "display_name": "Aorta",
    "body_region": "chest",
    "system": "cardiovascular",
    "snomed_code": "15825003",
    "fma_id": "FMA3734",
    "vista3d_label": 52,
    "structures": ["ascending aorta", "aortic arch", "descending aorta"],
    "common_pathologies": ["aneurysm", "dissection", "coarctation"],
    # ...
}
```

### 5.5 Keyword Routing

The `_get_knowledge_context()` method in `ImagingRAGEngine` scans the query text against all three dictionaries:

```python
def _get_knowledge_context(self, query: str) -> str:
    query_lower = query.lower()
    contexts = []

    for key in IMAGING_PATHOLOGIES:
        if key.replace("_", " ") in query_lower:
            contexts.append(get_pathology_context(key))
            break

    for key in IMAGING_MODALITIES:
        if key in query_lower:
            contexts.append(get_modality_context(key))
            break

    for key in IMAGING_ANATOMY:
        if key in query_lower:
            contexts.append(get_anatomy_context(key))
            break

    return "\n\n".join(contexts) if contexts else ""
```

The routing uses simple substring matching with underscore-to-space normalization. It breaks after the first match per domain, so only one pathology, one modality, and one anatomy entry are injected per query. This prevents knowledge context from overwhelming the LLM prompt.

### 5.6 Comparison Entity Resolution

The `resolve_comparison_entity()` and `get_comparison_context()` functions support the comparative analysis pipeline. They check all three dictionaries for entity matches and produce structured comparison tables.

---

## Chapter 6: Query Expansion Engineering

### 6.1 Architecture

Query expansion in `src/query_expansion.py` (186 lines) uses 12 domain-specific expansion maps:

| Map | Keys | Example Key -> Synonyms |
|---|---|---|
| `MODALITY_EXPANSION` | 9 | `"ct"` -> `{"computed tomography", "cat scan", "helical ct", ...}` |
| `BODY_REGION_EXPANSION` | 12 | `"chest"` -> `{"thorax", "thoracic", "lung", "pulmonary", ...}` |
| `PATHOLOGY_EXPANSION` | 22 | `"hemorrhage"` -> `{"bleeding", "hematoma", "ich", ...}` |
| `AI_TASK_EXPANSION` | 8 | `"segmentation"` -> `{"segment", "delineation", "contour", ...}` |
| `SEVERITY_EXPANSION` | 3 | `"critical"` -> `{"emergent", "stat", "life-threatening", ...}` |
| `FINDING_EXPANSION` | 6 | `"consolidation"` -> `{"airspace opacity", "air bronchograms", ...}` |
| `GUIDELINE_EXPANSION` | 9 | `"lung_rads"` -> `{"lung-rads", "ldct screening", ...}` |
| `DEVICE_EXPANSION` | 5 | `"510k"` -> `{"510(k)", "premarket notification", ...}` |
| `DATASET_EXPANSION` | 7 | `"lidc"` -> `{"lidc-idri", "lung image database", ...}` |
| `MODEL_ARCHITECTURE_EXPANSION` | 11 | `"vista3d"` -> `{"vista-3d", "versatile imaging segmentation", ...}` |
| `MEASUREMENT_EXPANSION` | 5 | `"hounsfield"` -> `{"hu", "density", "attenuation", ...}` |
| `CONTRAST_EXPANSION` | 4 | `"gadolinium"` -> `{"gad", "gd", "gbca", "mri contrast", ...}` |

### 6.2 The expand_query() Function

```python
def expand_query(query: str) -> Set[str]:
    query_lower = query.lower()
    expanded = set()

    for expansion_map in ALL_EXPANSION_MAPS:
        for key, synonyms in expansion_map.items():
            if key.replace("_", " ") in query_lower or key in query_lower:
                expanded.update(synonyms)

    return expanded
```

The function returns a set of expansion terms. The RAG engine appends up to 10 terms to the original query:

```python
expanded_terms = expand_query(query)
search_text = query
if expanded_terms:
    search_text = f"{query} {' '.join(list(expanded_terms)[:10])}"
```

### 6.3 Adding New Expansion Maps

To add a new expansion domain:

1. Define the map as a module-level dict:

```python
PROTOCOL_EXPANSION = {
    "low_dose": {"low dose ct", "ldct", "dose reduction", "low radiation"},
    "high_resolution": {"hrct", "thin slice", "1mm", "high resolution ct"},
}
```

2. Add it to `ALL_EXPANSION_MAPS`:

```python
ALL_EXPANSION_MAPS = [
    # ... existing maps ...
    PROTOCOL_EXPANSION,
]
```

### 6.4 Synonym Expansion vs Modality/Region Context

Query expansion and knowledge context serve different purposes:

- **Query expansion** modifies the embedding input to improve vector recall (same concept, different terminology)
- **Knowledge context** injects structured domain knowledge into the LLM prompt for synthesis accuracy

Both fire independently during `retrieve()`.

---

## Chapter 7: The Comparative Analysis System

### 7.1 Detection: _COMPARATIVE_RE

The regex pattern detects comparative intent:

```python
_COMPARATIVE_RE = re.compile(
    r"\b(compare|compared to|vs\.?|versus|difference between"
    r"|head.to.head|better than|advantages|disadvantages)\b",
    re.IGNORECASE,
)
```

Any query matching this pattern gets routed to the comparative pipeline instead of standard retrieval.

### 7.2 Entity Parsing: _ENTITY_PATTERNS

Eight regex patterns extract the two entities from a comparative query:

| Pattern | Example Match |
|---|---|
| `(.+?) vs\.? (.+)` | "CT vs MRI" |
| `(.+?) versus (.+)` | "CT versus MRI" |
| `compare (.+?) (?:and\|with\|to) (.+)` | "compare CT and MRI" |
| `(.+?) compared to (.+)` | "CT compared to MRI" |
| `difference(?:s)? between (.+?) and (.+)` | "difference between CT and MRI" |
| `(.+?) head[\s-]to[\s-]head (.+)` | "CT head-to-head MRI" |
| `advantages of (.+?) over (.+)` | "advantages of CT over MRI" |
| `(.+?) better than (.+)` | "CT better than MRI" |

Patterns are tried in order; the first match wins. A fallback splits on separator keywords if no regex matches.

### 7.3 Entity Resolution: IMAGING_ENTITY_MAP

The `IMAGING_ENTITY_MAP` contains 20+ entity definitions across five categories:

| Category | Entities |
|---|---|
| Modalities | ct, mri, pet, ultrasound, xray, mammography |
| Architectures | cnn, transformer, unet |
| Techniques | dlir, iterative |
| Tasks | detection, segmentation, classification |
| Technologies | photon_counting, dual_energy, conventional_ct |

Resolution uses a 5-step cascade:
1. Exact key match in `IMAGING_ENTITY_MAP`
2. Exact alias match (case-insensitive)
3. Fuzzy substring match in aliases
4. Fall back to `knowledge.py` entity resolution
5. Return raw text as `unknown` type

### 7.4 retrieve_comparative() Pipeline

```python
def retrieve_comparative(self, question, **kwargs) -> ComparativeResult:
    # 1. Parse entities
    entity_a_str, entity_b_str = self._parse_comparison_entities(question)

    # 2. Resolve entities against IMAGING_ENTITY_MAP
    entity_a_resolved = _resolve_imaging_entity(entity_a_str)
    entity_b_resolved = _resolve_imaging_entity(entity_b_str)

    # 3. Dual retrieval -- augment queries with full entity names
    search_a = f"{question} {entity_a_resolved.get('full_name', entity_a_str)}"
    search_b = f"{question} {entity_b_resolved.get('full_name', entity_b_str)}"
    evidence_a = self.retrieve(search_a, **kwargs)
    evidence_b = self.retrieve(search_b, **kwargs)

    # 4. Find shared evidence (same record_id in both result sets)
    shared = self._find_shared_evidence(evidence_a.hits, evidence_b.hits)

    # 5. Build knowledge-graph comparison context
    comparison_context = get_comparison_context(kg_a, kg_b)

    return ComparativeResult(...)
```

### 7.5 Dual-Query Evidence

The key insight is that each entity gets its own retrieval pass. The query for entity A is the original question augmented with entity A's full name, and similarly for entity B. This biases each search toward evidence relevant to that specific entity.

Shared evidence (records appearing in both result sets) is surfaced separately, as it often contains direct comparison studies.

### 7.6 Comparative System Prompt

The `COMPARATIVE_SYSTEM_PROMPT` instructs the LLM to produce a structured 7-section comparison:

1. Technical Specifications (table)
2. Clinical Performance (sensitivity, specificity, AUC)
3. Radiation Dose / Safety
4. AI Integration Capabilities
5. Clinical Guidelines & Recommendations
6. Cost-Effectiveness Considerations
7. Summary Recommendation

---

## Chapter 8: Clinical Workflow Development

### 8.1 Base Class Architecture

`src/workflows/base.py` (83 lines) defines `BaseImagingWorkflow`:

```python
class BaseImagingWorkflow(ABC):
    WORKFLOW_NAME: str = "base"
    TARGET_LATENCY_SEC: float = 60.0
    MODALITY: str = ""
    BODY_REGION: str = ""
    MODELS_USED: List[str] = []

    def __init__(self, mock_mode=True, nim_clients=None, mock_overrides=None):
        ...

    @abstractmethod
    def preprocess(self, input_path: str) -> Any: ...

    @abstractmethod
    def infer(self, preprocessed: Any) -> Dict: ...

    @abstractmethod
    def postprocess(self, inference_result: Dict) -> WorkflowResult: ...

    @abstractmethod
    def _mock_inference(self) -> Dict: ...

    def run(self, input_path: str = "") -> WorkflowResult: ...
```

### 8.2 The preprocess -> infer -> postprocess Pattern

**preprocess()** -- loads and prepares input data:
- CT workflows: load DICOM/NIfTI, apply windowing, resample to isotropic spacing
- CXR workflow: load DICOM, resize to model input dimensions, normalize pixel values
- MRI workflows: load multi-sequence data, apply bias field correction, register sequences

**infer()** -- runs model inference:
- Calls NIM clients (VISTA-3D, MAISI, VILA-M3) or MONAI models
- Returns raw inference dict with predictions, probabilities, segmentation masks

**postprocess()** -- extracts clinical meaning:
- Converts raw model outputs to findings, measurements, classifications
- Applies clinical scoring systems (Lung-RADS, BI-RADS, etc.)
- Determines severity level
- Returns a structured `WorkflowResult`

### 8.3 The run() Orchestrator

```python
def run(self, input_path: str = "") -> WorkflowResult:
    start = time.time()
    try:
        if self.mock_mode:
            raw = self._mock_inference()
            if self.mock_overrides:
                raw.update(self.mock_overrides)
        else:
            preprocessed = self.preprocess(input_path)
            raw = self.infer(preprocessed)

        result = self.postprocess(raw)
        result.inference_time_ms = (time.time() - start) * 1000
        result.is_mock = self.mock_mode
        result.workflow_name = self.WORKFLOW_NAME
        return result
    except Exception as e:
        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.FAILED,
            inference_time_ms=(time.time() - start) * 1000,
            is_mock=self.mock_mode,
        )
```

The `mock_overrides` parameter allows demo cases to inject specific findings (e.g., force Lung-RADS 4B classification for a demo scenario).

### 8.4 Six Implemented Workflows

| Workflow | File | Lines | Modality | Scoring System | Models |
|---|---|---|---|---|---|
| CT Head Hemorrhage | `ct_head_hemorrhage.py` | 488 | CT | Volume/midline shift | 3D U-Net, VISTA-3D |
| CT Lung Nodule | `ct_chest_lung_nodule.py` | 598 | CT | Lung-RADS v2022 | DenseNet, VISTA-3D |
| CT Coronary Angiography | `ct_coronary_angiography.py` | 357 | CT | CAD-RADS 2.0 | SegResNet |
| CXR Rapid Findings | `cxr_rapid_findings.py` | 658 | CXR | Multi-finding classification | DenseNet-121 |
| MRI Brain MS Lesion | `mri_brain_ms_lesion.py` | 645 | MRI | Lesion count/volume | nnU-Net, FLAIR |
| MRI Prostate PI-RADS | `mri_prostate_pirads.py` | 280 | MRI | PI-RADS v2.1 | SegResNet |

### 8.5 Adding a New Workflow

To add a PET/CT lymphoma staging workflow:

```python
# src/workflows/pet_ct_lymphoma_staging.py

from src.workflows.base import BaseImagingWorkflow
from src.models import FindingSeverity, WorkflowResult, WorkflowStatus

class PETCTLymphomaWorkflow(BaseImagingWorkflow):
    WORKFLOW_NAME = "pet_ct_lymphoma_staging"
    TARGET_LATENCY_SEC = 120.0
    MODALITY = "pet_ct"
    BODY_REGION = "whole_body"
    MODELS_USED = ["VISTA-3D", "SUV quantification"]

    def preprocess(self, input_path: str):
        # Load PET/CT DICOM series
        # Register PET to CT
        # Extract SUV maps
        return {"pet_volume": pet_vol, "ct_volume": ct_vol}

    def infer(self, preprocessed):
        # Run organ segmentation via VISTA-3D
        # Identify FDG-avid regions
        # Calculate SUVmax per region
        return {"regions": [...], "suv_values": {...}}

    def postprocess(self, inference_result) -> WorkflowResult:
        # Apply Deauville scoring (1-5)
        # Map to Lugano classification
        # Determine staging
        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            classification=f"Deauville {score}",
            severity=self._deauville_to_severity(score),
            findings=findings,
            measurements={"suv_max": max_suv},
        )

    def _mock_inference(self):
        return {
            "regions": [{"name": "mediastinal", "suv_max": 8.2}],
            "deauville_score": 4,
        }
```

Register it in `src/workflows/__init__.py`:

```python
from .pet_ct_lymphoma_staging import PETCTLymphomaWorkflow

WORKFLOW_REGISTRY["pet_ct_lymphoma_staging"] = PETCTLymphomaWorkflow
```

Add routing rules in `api/routes/events.py`:

```python
WORKFLOW_ROUTING[("PT", "whole_body")] = "pet_ct_lymphoma_staging"
```

---

## Chapter 9: NIM Integration Deep Dive

### 9.1 BaseNIMClient Architecture

`src/nim/base.py` (83 lines) provides the foundation for all NIM clients:

```python
class BaseNIMClient(ABC):
    def __init__(self, base_url, service_name, mock_enabled=True):
        self.base_url = base_url.rstrip("/")
        self.service_name = service_name
        self.mock_enabled = mock_enabled
        self._available = None
        self._last_check = 0
        self._check_interval = 30.0  # Cached health check interval
```

Key methods:
- `health_check()` -- HTTP GET to `/v1/health/ready`
- `is_available()` -- cached health check (30-second TTL)
- `_request()` -- HTTP POST with tenacity retry (3 attempts, exponential backoff)
- `_invoke_or_mock()` -- try real NIM, fall back to mock if unavailable

### 9.2 The Fallback Chain

Each NIM client implements a multi-tier fallback:

**LLM Client (`llm_client.py`, 439 lines):**
```
Local NIM (Llama-3-70B) -> NVIDIA Cloud NIM (Llama-3.1-8B) -> Anthropic Claude -> Mock
```

**VILA-M3 Client (`vilam3_client.py`, 537 lines):**
```
Local VILA-M3 NIM -> NVIDIA Cloud NIM (Llama-3.2-11B-Vision) -> Mock
```

**VISTA-3D Client (`vista3d_client.py`, 215 lines):**
```
Local VISTA-3D NIM -> Mock (realistic segmentation results)
```

**MAISI Client (`maisi_client.py`, 136 lines):**
```
Local MAISI NIM -> Mock (synthetic CT metadata)
```

### 9.3 Client Architecture Details

**LlamaLLMClient** uses the OpenAI Python client because NIM exposes an OpenAI-compatible `/v1/chat/completions` endpoint:

```python
class LlamaLLMClient(BaseNIMClient):
    def __init__(self, base_url, mock_enabled=True,
                 anthropic_api_key=None, nvidia_api_key=None,
                 cloud_url="https://integrate.api.nvidia.com/v1",
                 cloud_llm_model="meta/llama-3.1-8b-instruct",
                 local_llm_model="meta/llama3-70b-instruct"):
```

The health check for the LLM client uses `/v1/models` instead of `/v1/health/ready` because the OpenAI-compatible API exposes model listing.

**VILAM3Client** supports multimodal inputs by base64-encoding images and sending them in the OpenAI vision message format:

```python
def _encode_image(self, image_path: str) -> str:
    path = Path(image_path)
    return base64.b64encode(path.read_bytes()).decode("utf-8")
```

### 9.4 Mock Mode

Mock mode is critical for development, testing, and demos without GPU access. Each client implements `_mock_response()` returning clinically realistic synthetic data:

- VISTA-3D mock: 3 segmented organs with realistic volumes
- MAISI mock: 512x512x512 volume metadata, 127 annotated classes
- VILA-M3 mock: structured radiology report text
- LLM mock: generic clinical assessment text

The mock system is controlled by:
- `NIM_ALLOW_MOCK_FALLBACK: bool = True` in settings
- Per-client `mock_enabled` constructor parameter
- `NIM_MODE: str = "local"` (values: "local", "cloud", "mock")

### 9.5 NIMServiceManager

The `NIMServiceManager` aggregates all four NIM clients and provides:
- `check_all_services()` -- returns status dict for all services
- `get_available_services()` -- list of services with live endpoints
- `get_client(name)` -- retrieve a specific client by name

---

## Chapter 10: Testing Strategies

### 10.1 Test Suite Overview

The test suite consists of 12 modules with 620 tests, all passing:

| Module | Tests | Focus |
|---|---|---|
| `test_models.py` | ~80 | Pydantic model validation, enum coverage |
| `test_export.py` | ~60 | Markdown, JSON, PDF, FHIR R4 output |
| `test_rag_engine.py` | ~75 | Retrieve, query, comparative, scoring |
| `test_workflows.py` | ~65 | All 6 workflows, mock inference, postprocess |
| `test_agent.py` | ~55 | Agent orchestrator, search planning |
| `test_knowledge.py` | ~50 | Knowledge graph lookups, comparison context |
| `test_nim_clients.py` | ~60 | All 4 NIM clients, fallback chains |
| `test_dicom_ingestion.py` | ~40 | DICOM watcher, Orthanc events |
| `test_cross_modal.py` | ~45 | Trigger evaluation, genomic queries |
| `test_query_expansion.py` | ~35 | All 12 expansion maps |
| `test_preview_generator.py` | ~30 | Slice animation, windowing |
| `test_preview_api.py` | ~25 | Preview API endpoints |

### 10.2 The Mock-Everything Approach

Every external dependency is mocked in `tests/conftest.py`. The test suite runs without:
- Milvus database
- NIM GPU services
- Network access (PubMed, ClinicalTrials.gov)
- Orthanc DICOM server
- Anthropic/NVIDIA API keys

### 10.3 Core Fixtures

`conftest.py` provides 10+ fixtures:

**mock_embedder** -- returns 384-dim numpy vectors (random):

```python
@pytest.fixture
def mock_embedder():
    embedder = MagicMock()

    def _encode(texts, normalize_embeddings=True, **kwargs):
        if isinstance(texts, str):
            return np.random.randn(384).astype(np.float32)
        return np.random.randn(len(texts), 384).astype(np.float32)

    embedder.encode = MagicMock(side_effect=_encode)
    return embedder
```

**mock_llm_client** -- returns template clinical text:

```python
@pytest.fixture
def mock_llm_client():
    client = MagicMock()
    client.generate.return_value = (
        "Based on the available imaging evidence, the findings are consistent "
        "with normal anatomy without acute pathology."
    )
    client.generate_stream.return_value = iter([
        "Based on ", "the available ", "imaging evidence, ", ...
    ])
    return client
```

**mock_collection_manager** -- returns 3 sample search results per collection:

```python
@pytest.fixture
def mock_collection_manager():
    manager = MagicMock()
    manager.search.return_value = [
        {"id": "lit-001", "score": 0.85, "text_chunk": "AI-based hemorrhage detection..."},
        {"id": "lit-002", "score": 0.78, "text_chunk": "VISTA-3D segments 127 structures..."},
        {"id": "lit-003", "score": 0.72, "text_chunk": "Lung-RADS v2022 classifies..."},
    ]
    manager.get_collection_stats.return_value = {name: 100 for name in IMAGING_COLLECTION_NAMES}
    return manager
```

**mock_nim_services** -- all 4 NIM clients with mock responses (VISTA-3D segmentation, MAISI generation, VILA-M3 image analysis, LLM generation).

**sample_search_hits** -- 5 SearchHit objects spanning literature, trials, findings, guidelines, and benchmarks.

**sample_evidence** -- a `CrossCollectionResult` with 5 hits and knowledge context.

**sample_workflow_result** -- a `WorkflowResult` from CT head hemorrhage (urgent, 12.5 mL volume).

**sample_agent_response** -- complete `AgentResponse` for export testing.

### 10.4 Testing Patterns

**Workflow tests** verify postprocess logic by providing mock inference results and checking classifications:

```python
def test_lung_nodule_4a(mock_lung_nodule_workflow):
    result = mock_lung_nodule_workflow.run()
    assert "4A" in result.classification or "Lung-RADS" in result.classification
    assert result.status == WorkflowStatus.COMPLETED
```

**RAG engine tests** verify weighted scoring:

```python
def test_retrieve_weighted_scoring(rag_engine):
    result = rag_engine.retrieve("hemorrhage detection CT")
    # Literature weight (0.18) > Findings weight (0.15)
    # So literature hits should rank higher for equal raw scores
    assert result.hits[0].collection == "imaging_literature"
```

**Cross-modal tests** verify trigger thresholds:

```python
def test_lung_rads_4b_triggers(trigger, mock_embedder):
    result = WorkflowResult(
        workflow_name="ct_chest_lung_nodule",
        classification="Lung-RADS 4B",
        severity=FindingSeverity.URGENT,
    )
    cross_modal = trigger.evaluate(result)
    assert cross_modal is not None
    assert "Lung-RADS 4B" in cross_modal.trigger_reason
```

---

## Chapter 11: Export System Deep Dive

### 11.1 Four Export Formats

`src/export.py` (691 lines) provides four export functions:

| Function | Output | Use Case |
|---|---|---|
| `export_markdown()` | String | Chat display, clipboard, documentation |
| `export_json()` | String | API responses, downstream integration |
| `export_pdf()` | File path | Clinical reports, presentations |
| `export_fhir_r4()` | Dict | EHR integration, interoperability |

### 11.2 Markdown Export

```python
def export_markdown(response: AgentResponse) -> str:
```

Produces a structured Markdown report with sections:
- Header with query and timestamp
- Analysis (the LLM-synthesized answer)
- Evidence table (grouped by collection, top 5 per collection)
- Workflow Results (findings, measurements, classification)
- NIM Services Used
- Research-use-only disclaimer

### 11.3 PDF Export with ReportLab

```python
def export_pdf(response: AgentResponse, output_path: str) -> str:
```

Generates NVIDIA-branded PDF reports using ReportLab:
- **Header:** NVIDIA green (RGB 118, 185, 0) bar with title
- **Page setup:** Letter size, 40pt margins
- **Sections:** Query, Analysis (markdown-to-PDF conversion), Evidence table, Workflow results
- **Severity colors:** Critical=red, Urgent=orange, Significant=amber, Routine=green, Normal=gray

The `_clean_markdown()` helper strips Markdown formatting for PDF paragraph rendering, converting `**bold**` to `<b>bold</b>` tags compatible with ReportLab's `Paragraph`.

### 11.4 FHIR R4 Export

The FHIR R4 export produces a standards-compliant `DiagnosticReport` resource with embedded `Observation` resources:

**Code systems used:**

| System | URI | Purpose |
|---|---|---|
| LOINC | `http://loinc.org` | Observation codes |
| SNOMED CT | `http://snomed.info/sct` | Finding codes |
| DICOM | `http://dicom.nema.org/resources/ontology/DCM` | Modality codes |
| HL7 Interpretation | `http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation` | Severity mapping |

**SNOMED finding codes:**

| Finding | SNOMED Code |
|---|---|
| hemorrhage | 50960005 |
| nodule | 416940007 |
| consolidation | 95436008 |
| pneumothorax | 36118008 |
| effusion | 60046008 |
| fracture | 125605004 |
| cardiomegaly | 8186001 |
| mass | 4147007 |
| edema | 267038008 |
| normal | 17621005 |

**Severity to FHIR Interpretation mapping:**

| Severity | FHIR Code |
|---|---|
| critical | HH (critically high) |
| urgent | H (high) |
| significant | A (abnormal) |
| routine | N (normal) |
| normal | N (normal) |

**DICOM modality codes:**

| Agent Modality | DICOM Code |
|---|---|
| ct | CT |
| mri | MR |
| xray | DX |
| cxr | CR |
| ultrasound | US |
| pet | PT |
| mammography | MG |
| fluoroscopy | RF |

---

## Chapter 12: Cross-Modal Integration

### 12.1 Architecture

The cross-modal system in `src/cross_modal.py` (401 lines) bridges imaging findings to the genomic_evidence collection (3.5M vectors) populated by the rag-chat-pipeline.

```python
class CrossModalTrigger:
    GENOMIC_COLLECTION = "genomic_evidence"
    TOP_K_PER_QUERY = 5
    SCORE_THRESHOLD = 0.40
```

### 12.2 Five Evaluator Functions

Each evaluator maps a workflow result to genomic query templates:

| Evaluator | Workflow | Trigger Condition | Genomic Queries |
|---|---|---|---|
| `_evaluate_lung_nodule` | `ct_chest_lung_nodule` | Lung-RADS 4A, 4B, 4X | EGFR, ALK, ROS1, KRAS driver mutations |
| `_evaluate_cxr` | `cxr_rapid_findings` | Critical/Urgent + consolidation | Antimicrobial resistance, immune response |
| `_evaluate_brain_lesion` | `mri_brain_ms_lesion` | Classification contains "highly_active" | HLA-DRB1 susceptibility, treatment response |
| `_evaluate_cardiac` | `ct_coronary_angiography` | CAD-RADS >= 3 | LDLR, PCSK9, APOB, polygenic risk |
| `_evaluate_prostate` | `mri_prostate_pirads` | PI-RADS >= 4 | BRCA2, HOXB13, ATM susceptibility |

### 12.3 Severity Threshold Mapping

```
Workflow Result -> Evaluator -> Severity Check -> Genomic Queries -> Milvus Search -> CrossModalResult
```

Each evaluator parses the `classification` string from the `WorkflowResult` and checks against specific high-risk categories:

```python
# Lung nodule: parse Lung-RADS category
high_risk_categories = ["4A", "4B", "4X"]
high_risk = any(cat in classification for cat in high_risk_categories)

# Cardiac CTA: parse CAD-RADS score
high_risk_categories = ["3", "4A", "4B", "5"]
high_risk = any(f"CAD-RADS {cat}" in classification for cat in high_risk_categories)

# Prostate MRI: parse PI-RADS score
high_risk_scores = ["4", "5"]
high_risk = any(f"PI-RADS {score}" in classification for score in high_risk_scores)
```

### 12.4 Genomic Query Templates

Each trigger fires 2-3 pre-defined genomic queries:

```python
LUNG_CANCER_QUERIES = [
    "lung cancer driver mutations EGFR ALK ROS1 KRAS",
    "non-small cell lung cancer NSCLC targeted therapy genomics",
    "lung adenocarcinoma molecular subtypes precision medicine",
]

CARDIAC_GENOMICS_QUERIES = [
    "cardiovascular disease genetic risk factors LDLR PCSK9 APOB",
    "familial hypercholesterolemia genomic variants",
    "coronary artery disease polygenic risk score",
]
```

### 12.5 The _query_genomics() Method

```python
def _query_genomics(self, queries, trigger_reason):
    for query in queries:
        embedding = self.embedder.encode(query, normalize_embeddings=True)
        hits = self.collection_manager.search(
            collection_name=self.GENOMIC_COLLECTION,
            query_embedding=embedding_list,
            top_k=self.TOP_K_PER_QUERY,      # 5 per query
            score_threshold=self.SCORE_THRESHOLD,  # 0.40 minimum
        )
        # De-duplicate by ID across queries
```

Results are de-duplicated by record ID across all queries and assembled into a `CrossModalResult` with genomic context strings, hit count, and an enrichment summary.

### 12.6 Integration with Agent

The `ImagingIntelligenceAgent` invokes cross-modal triggers after workflow execution:

```python
if self.cross_modal_trigger and workflow_result:
    cross_modal_result = self.cross_modal_trigger.evaluate(workflow_result)
    if cross_modal_result:
        response.cross_modal = cross_modal_result
```

The cross-modal enrichment is controlled by the `CROSS_MODAL_ENABLED` setting (default: `False` for safety in production).

---

## Chapter 13: Production Deployment

### 13.1 Docker Multi-Stage Build

The `Dockerfile` uses a multi-stage build:

**Stage 1: Builder** -- installs Python dependencies and compiles native extensions

**Stage 2: Runtime** -- copies only the installed packages and application code

### 13.2 Compose Topology

**Full stack** (`docker-compose.yml`) -- 11 services:

```
orthanc (8042, 4242)          -- DICOM server + web viewer
ohif-viewer (8526)            -- OHIF DICOM viewer
milvus-etcd                   -- Milvus metadata (etcd)
milvus-minio                  -- Milvus object storage (MinIO)
milvus-standalone (19530)     -- Milvus vector database
imaging-streamlit (8525)      -- Streamlit chat UI
imaging-api (8524)            -- FastAPI REST server
imaging-setup                 -- One-shot collection + seed
nim-llm (8520)                -- Meta Llama 3 NIM
nim-vista3d (8530)            -- NVIDIA VISTA-3D NIM
nim-maisi (8531)              -- NVIDIA MAISI NIM
nim-vilam3 (8532)             -- VILA-M3 VLM NIM
```

**Lite stack** (`docker-compose.lite.yml`) -- 6 services (no GPU required):

```
milvus-etcd                   -- Milvus metadata
milvus-minio                  -- Milvus object storage
milvus-standalone (19530)     -- Milvus vector database
imaging-streamlit (8525)      -- Streamlit chat UI
imaging-api (8524)            -- FastAPI REST server
imaging-setup                 -- One-shot collection + seed
```

All NIM-dependent features run in mock mode with the Lite stack.

### 13.3 Port Map

| Service | Port | Protocol |
|---|---|---|
| Streamlit UI | 8525 | HTTP |
| FastAPI REST | 8524 | HTTP |
| Milvus gRPC | 19530 | gRPC |
| Milvus metrics | 9091 | HTTP |
| Orthanc REST/Web | 8042 | HTTP |
| Orthanc DICOM | 4242 | DICOM C-STORE |
| OHIF Viewer | 8526 | HTTP |
| NIM LLM | 8520 | HTTP (OpenAI-compat) |
| NIM VISTA-3D | 8530 | HTTP |
| NIM MAISI | 8531 | HTTP |
| NIM VILA-M3 | 8532 | HTTP |

### 13.4 Health Checks

**Orthanc:**
```yaml
healthcheck:
  test: ["CMD", "curl", "-f", "http://localhost:8042/system"]
  interval: 30s
  timeout: 10s
  retries: 5
```

**Milvus:**
```yaml
healthcheck:
  test: ["CMD", "curl", "-f", "http://localhost:9091/healthz"]
  interval: 30s
  timeout: 20s
  retries: 3
```

**FastAPI:**
The `/health` endpoint returns service status and Milvus connectivity.

**NIM services:**
Each NIM client checks `/v1/health/ready` (or `/v1/models` for LLM) with a 30-second cached TTL.

### 13.5 Monitoring

**Prometheus metrics** exposed by the FastAPI server:

| Metric | Type | Labels | Description |
|---|---|---|---|
| `imaging_agent_queries_total` | Counter | `endpoint` | Total RAG queries |
| `imaging_agent_query_duration_seconds` | Histogram | `endpoint` | Query latency |
| `imaging_agent_search_hits` | Histogram | -- | Evidence hits per query |

Histogram buckets for latency: 0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0 seconds.

### 13.6 Environment Variables

All settings use the `IMAGING_` prefix:

```bash
IMAGING_MILVUS_HOST=milvus-standalone
IMAGING_MILVUS_PORT=19530
IMAGING_ANTHROPIC_API_KEY=sk-ant-...
IMAGING_NIM_LLM_URL=http://nim-llm:8520/v1
IMAGING_NIM_MODE=local
IMAGING_CROSS_MODAL_ENABLED=true
```

### 13.7 DICOM Event Integration

The `api/routes/events.py` module provides webhook endpoints for Orthanc. When a DICOM study arrives:

1. Orthanc fires a change event to the webhook
2. The event handler extracts modality and body region from DICOM tags
3. `WORKFLOW_ROUTING` maps `(modality, body_region)` to a workflow name:

```python
WORKFLOW_ROUTING = {
    ("CT", "head"): "ct_head_hemorrhage",
    ("CT", "brain"): "ct_head_hemorrhage",
    ("CT", "chest"): "ct_chest_lung_nodule",
    ("CT", "lung"): "ct_chest_lung_nodule",
    ("CR", "chest"): "cxr_rapid_findings",
    ("DX", "chest"): "cxr_rapid_findings",
    ("MR", "brain"): "mri_brain_ms_lesion",
    ("MR", "head"): "mri_brain_ms_lesion",
}
```

4. The matched workflow runs and results are stored in an in-memory deque (max 200 entries)

---

## Chapter 14: Future Architecture

### 14.1 Federated Learning with NVIDIA FLARE

The `flare/` directory contains federated learning job configurations for three clinical tasks:

| Job | Directory | Description |
|---|---|---|
| CXR Classification | `flare/job_configs/cxr_classification/` | Multi-institution chest X-ray classification |
| CT Segmentation | `flare/job_configs/ct_segmentation/` | Federated organ segmentation |
| Lung Nodule Detection | `flare/job_configs/lung_nodule_detection/` | Distributed lung nodule CADe |

FLARE enables training across hospital sites without centralizing patient data. Each site trains locally and shares only model weight updates (federated averaging).

**Integration pattern:**
```
Site 1 (Hospital A)     Site 2 (Hospital B)     Site 3 (Hospital C)
     |                       |                       |
     v                       v                       v
Local Training          Local Training          Local Training
     |                       |                       |
     +--------> FLARE Server <--------+
                    |
                    v
            Aggregated Model
                    |
                    v
           Deploy to NIM Endpoint
```

### 14.2 Multi-Agent Systems

The Imaging Intelligence Agent is designed to operate within the HCLS AI Factory's multi-agent ecosystem:

- **Biomarker Discovery Agent** -- variant interpretation, biomarker identification
- **Precision Oncology Agent** -- treatment recommendation, clinical trial matching
- **CAR-T Intelligence Agent** -- cell therapy optimization
- **Autoimmune Intelligence Agent** -- autoimmune disease analysis

Future cross-agent integration patterns:
- Imaging findings triggering oncology treatment recommendations
- Genomic variants informing imaging protocol selection
- Federated knowledge sharing across agent collections

### 14.3 Fine-Tuned Embeddings

The current BGE-small-en-v1.5 model is a general-purpose text encoder. Future improvements include:

**Domain-adapted embeddings:**
- Fine-tune on RadLex terminology and radiology reports
- Use contrastive learning on (query, relevant_passage) pairs from radiology literature
- Target improvement: +5-10% recall on radiology-specific queries

**Multi-modal embeddings:**
- Encode DICOM image features alongside text
- Use BiomedCLIP or RadFM for joint image-text embeddings
- Enable "find studies similar to this image" searches

### 14.4 Scalability Roadmap

| Current | Near-Term | Long-Term |
|---|---|---|
| IVF_FLAT, 11 collections | HNSW, 15+ collections | Distributed Milvus, 50+ collections |
| BGE-small (384d) | Fine-tuned BGE (384d) | Multi-modal (768d+) |
| Mock NIM fallback | Full NIM deployment | Multi-GPU NIM cluster |
| Single DGX Spark | DGX Spark + Cloud NIM | Multi-node DGX |
| In-memory event history | PostgreSQL event store | Kafka event streaming |

---

## Appendix A: Complete API Reference

### A.1 Core Query Endpoints

**POST /query**

RAG query with multi-collection search, knowledge augmentation, and LLM synthesis.

Request:
```json
{
    "question": "What is the sensitivity of AI for hemorrhage detection on CT?",
    "modality": "ct",
    "body_region": "head",
    "top_k": 5,
    "include_genomic": true,
    "include_nim": true,
    "collections": null,
    "year_min": 2020,
    "year_max": null,
    "conversation_context": ""
}
```

Response:
```json
{
    "question": "...",
    "answer": "Based on the available evidence...",
    "evidence_count": 33,
    "collections_searched": 11,
    "search_time_ms": 45.2
}
```

**POST /query/stream**

Streaming version of `/query`. Returns Server-Sent Events (SSE) with incremental tokens.

**POST /query/comparative**

Explicitly invokes the comparative analysis pipeline for A-vs-B queries.

### A.2 Search Endpoints

**POST /search**

Raw vector search without LLM synthesis. Returns ranked evidence hits.

**POST /search/collection/{collection_name}**

Search a specific collection by name.

**GET /collections**

List all registered collections with record counts and status.

**GET /collections/{name}/stats**

Detailed statistics for a single collection.

### A.3 Workflow Endpoints

**POST /workflows/run**

Execute a named imaging workflow.

Request:
```json
{
    "workflow_name": "ct_head_hemorrhage",
    "input_path": "/data/studies/CT001",
    "mock_mode": true,
    "mock_overrides": {"classification": "urgent_hemorrhage"}
}
```

Response:
```json
{
    "workflow_name": "ct_head_hemorrhage",
    "status": "completed",
    "findings": [...],
    "measurements": {"volume_ml": 12.5, "midline_shift_mm": 3.2},
    "classification": "urgent_hemorrhage",
    "severity": "urgent",
    "inference_time_ms": 1250.0,
    "is_mock": true
}
```

**GET /workflows**

List available workflows with metadata.

**GET /workflows/{name}/info**

Workflow details: modality, body region, models used, target latency.

### A.4 NIM Service Endpoints

**GET /nim/status**

Status of all 4 NIM services (available/mock/unavailable).

Response:
```json
{
    "vista3d": "mock",
    "maisi": "mock",
    "vila_m3": "mock",
    "llm": "available"
}
```

**POST /nim/vista3d/segment**

Invoke VISTA-3D organ segmentation.

**POST /nim/maisi/generate**

Invoke MAISI synthetic CT generation.

**POST /nim/vilam3/analyze**

Invoke VILA-M3 image analysis with visual question answering.

### A.5 Report/Export Endpoints

**POST /reports/markdown**

Export an agent response as Markdown.

**POST /reports/json**

Export as structured JSON.

**POST /reports/pdf**

Export as NVIDIA-branded PDF file.

**POST /reports/fhir**

Export as FHIR R4 DiagnosticReport resource.

### A.6 Preview Endpoints

**POST /preview/generate**

Generate MP4/GIF preview from a DICOM series or NIfTI volume.

Request:
```json
{
    "study_id": "orthanc-study-001",
    "format": "mp4",
    "plane": "axial",
    "window": "lung",
    "fps": 8,
    "max_frames": 200
}
```

**GET /preview/{preview_id}**

Retrieve a previously generated preview file.

### A.7 Protocol/Dose Endpoints

**POST /protocol/recommend**

AI-optimized protocol recommendation for a clinical scenario.

**POST /dose/compare**

Compare standard vs AI-optimized radiation dose for a protocol.

### A.8 Demo Case Endpoints

**GET /demo/cases**

List available demo cases with metadata.

**POST /demo/run/{case_id}**

Execute a pre-built demo scenario end-to-end.

### A.9 Event Endpoints

**POST /events/dicom/study-complete**

Webhook for Orthanc study-complete events. Triggers automatic workflow routing.

**GET /events/history**

Retrieve recent DICOM ingestion history (max 200 entries).

### A.10 System Endpoints

**GET /health**

Health check with Milvus connectivity status.

**GET /metrics**

Prometheus metrics endpoint (query counts, latencies, hit distributions).

**GET /config**

Current configuration (sensitive values redacted).

---

## Appendix B: Configuration Reference

All environment variables use the `IMAGING_` prefix (controlled by Pydantic `env_prefix`).

### B.1 Path Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_PROJECT_ROOT` | Auto-detected | Project root directory |
| `IMAGING_DATA_DIR` | `{PROJECT_ROOT}/data` | Data storage directory |
| `IMAGING_CACHE_DIR` | `{DATA_DIR}/cache` | Cache directory |
| `IMAGING_REFERENCE_DIR` | `{DATA_DIR}/reference` | Reference data directory |
| `IMAGING_RAG_PIPELINE_ROOT` | `/home/adam/.../rag-chat-pipeline` | RAG pipeline root |

### B.2 Milvus Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_MILVUS_HOST` | `localhost` | Milvus server hostname |
| `IMAGING_MILVUS_PORT` | `19530` | Milvus gRPC port |
| `IMAGING_COLLECTION_LITERATURE` | `imaging_literature` | Literature collection name |
| `IMAGING_COLLECTION_TRIALS` | `imaging_trials` | Trials collection name |
| `IMAGING_COLLECTION_FINDINGS` | `imaging_findings` | Findings collection name |
| `IMAGING_COLLECTION_PROTOCOLS` | `imaging_protocols` | Protocols collection name |
| `IMAGING_COLLECTION_DEVICES` | `imaging_devices` | Devices collection name |
| `IMAGING_COLLECTION_ANATOMY` | `imaging_anatomy` | Anatomy collection name |
| `IMAGING_COLLECTION_BENCHMARKS` | `imaging_benchmarks` | Benchmarks collection name |
| `IMAGING_COLLECTION_GUIDELINES` | `imaging_guidelines` | Guidelines collection name |
| `IMAGING_COLLECTION_REPORT_TEMPLATES` | `imaging_report_templates` | Report templates collection |
| `IMAGING_COLLECTION_DATASETS` | `imaging_datasets` | Datasets collection name |
| `IMAGING_COLLECTION_GENOMIC` | `genomic_evidence` | Cross-agent genomic collection |

### B.3 Embedding Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | Embedding model name |
| `IMAGING_EMBEDDING_DIMENSION` | `384` | Vector dimension |
| `IMAGING_EMBEDDING_BATCH_SIZE` | `32` | Batch size for encoding |

### B.4 LLM Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_LLM_PROVIDER` | `anthropic` | LLM provider |
| `IMAGING_LLM_MODEL` | `claude-sonnet-4-6` | LLM model name |
| `IMAGING_ANTHROPIC_API_KEY` | `None` | Anthropic API key |

### B.5 NIM Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_NIM_LLM_URL` | `http://localhost:8520/v1` | Local LLM NIM endpoint |
| `IMAGING_NIM_VISTA3D_URL` | `http://localhost:8530` | VISTA-3D NIM endpoint |
| `IMAGING_NIM_MAISI_URL` | `http://localhost:8531` | MAISI NIM endpoint |
| `IMAGING_NIM_VILAM3_URL` | `http://localhost:8532` | VILA-M3 NIM endpoint |
| `IMAGING_NIM_MODE` | `local` | NIM mode: local, cloud, mock |
| `IMAGING_NIM_ALLOW_MOCK_FALLBACK` | `True` | Allow mock fallback |
| `IMAGING_NGC_API_KEY` | `None` | NGC API key for NIM pulls |
| `IMAGING_NIM_LOCAL_LLM_MODEL` | `meta/llama3-70b-instruct` | Local NIM LLM model |

### B.6 NVIDIA Cloud NIM Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_NVIDIA_API_KEY` | `None` | NVIDIA cloud API key |
| `IMAGING_NIM_CLOUD_URL` | `https://integrate.api.nvidia.com/v1` | Cloud NIM base URL |
| `IMAGING_NIM_CLOUD_LLM_MODEL` | `meta/llama-3.1-8b-instruct` | Cloud LLM model |
| `IMAGING_NIM_CLOUD_VLM_MODEL` | `meta/llama-3.2-11b-vision-instruct` | Cloud VLM model |

### B.7 RAG Search Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_TOP_K_PER_COLLECTION` | `5` | Results per collection |
| `IMAGING_SCORE_THRESHOLD` | `0.40` | Minimum similarity score |
| `IMAGING_WEIGHT_LITERATURE` | `0.18` | Literature collection weight |
| `IMAGING_WEIGHT_TRIALS` | `0.12` | Trials collection weight |
| `IMAGING_WEIGHT_FINDINGS` | `0.15` | Findings collection weight |
| `IMAGING_WEIGHT_PROTOCOLS` | `0.08` | Protocols collection weight |
| `IMAGING_WEIGHT_DEVICES` | `0.08` | Devices collection weight |
| `IMAGING_WEIGHT_ANATOMY` | `0.06` | Anatomy collection weight |
| `IMAGING_WEIGHT_BENCHMARKS` | `0.08` | Benchmarks collection weight |
| `IMAGING_WEIGHT_GUIDELINES` | `0.10` | Guidelines collection weight |
| `IMAGING_WEIGHT_REPORT_TEMPLATES` | `0.05` | Report templates weight |
| `IMAGING_WEIGHT_DATASETS` | `0.06` | Datasets collection weight |
| `IMAGING_WEIGHT_GENOMIC` | `0.04` | Genomic collection weight |

### B.8 Citation Scoring

| Variable | Default | Description |
|---|---|---|
| `IMAGING_CITATION_HIGH_THRESHOLD` | `0.75` | High-confidence citation threshold |
| `IMAGING_CITATION_MEDIUM_THRESHOLD` | `0.60` | Medium-confidence citation threshold |

### B.9 API Server Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_API_HOST` | `0.0.0.0` | API bind host |
| `IMAGING_API_PORT` | `8524` | API bind port |
| `IMAGING_STREAMLIT_PORT` | `8525` | Streamlit UI port |
| `IMAGING_API_BASE_URL` | `http://localhost:8524` | Internal API URL |
| `IMAGING_CORS_ORIGINS` | `http://localhost:8080,...` | CORS allowed origins |
| `IMAGING_MAX_REQUEST_SIZE_MB` | `10` | Max request body size |

### B.10 Orthanc/DICOM Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_ORTHANC_URL` | `http://localhost:8042` | Orthanc REST API URL |
| `IMAGING_ORTHANC_USERNAME` | `admin` | Orthanc username |
| `IMAGING_ORTHANC_PASSWORD` | `""` | Orthanc password |
| `IMAGING_DICOM_AUTO_INGEST` | `False` | Auto-ingest new DICOM studies |
| `IMAGING_DICOM_WATCH_INTERVAL` | `5` | Orthanc poll interval (seconds) |
| `IMAGING_DICOM_SERVER_URL` | `http://localhost:8042` | Legacy alias |

### B.11 OHIF Viewer Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_OHIF_URL` | `http://localhost:8526` | OHIF Viewer URL |

### B.12 Preview Generation Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_PREVIEW_CACHE_DIR` | `{DATA_DIR}/cache/previews` | Preview cache directory |
| `IMAGING_PREVIEW_DEFAULT_FPS` | `8` | Default frames per second |
| `IMAGING_PREVIEW_DEFAULT_FORMAT` | `mp4` | Default format (mp4 or gif) |
| `IMAGING_PREVIEW_MAX_FRAMES` | `200` | Maximum frames per preview |

### B.13 Data Source Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_NCBI_API_KEY` | `None` | PubMed API key (optional) |
| `IMAGING_PUBMED_MAX_RESULTS` | `5000` | Max PubMed results per query |
| `IMAGING_CT_GOV_BASE_URL` | `https://clinicaltrials.gov/api/v2` | ClinicalTrials.gov API |

### B.14 Scheduler/Ingest Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_INGEST_SCHEDULE_HOURS` | `168` | Ingest interval (168h = weekly) |
| `IMAGING_INGEST_ENABLED` | `False` | Enable scheduled ingestion |

### B.15 Cross-Modal Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_CROSS_MODAL_ENABLED` | `False` | Enable imaging->genomics triggers |

### B.16 Conversation Memory

| Variable | Default | Description |
|---|---|---|
| `IMAGING_MAX_CONVERSATION_CONTEXT` | `3` | Prior exchanges to inject |

### B.17 Prometheus Configuration

| Variable | Default | Description |
|---|---|---|
| `IMAGING_METRICS_ENABLED` | `True` | Enable Prometheus metrics |

---

## Appendix C: Scoring System Reference

### C.1 Lung-RADS (ACR Lung-RADS v2022)

Used by: `ct_chest_lung_nodule` workflow

| Category | Description | Finding Criteria | Management |
|---|---|---|---|
| 0 | Incomplete | Additional imaging needed | Complete evaluation |
| 1 | Negative | No nodules | Continue annual screening |
| 2 | Benign | Perifissural, calcified, or stable nodules | Continue annual screening |
| 3 | Probably benign | Solid nodule 6-8mm, part-solid <6mm | 6-month follow-up |
| 4A | Suspicious | Solid nodule 8-15mm, new part-solid 6mm+ | 3-month LDCT or PET/CT |
| 4B | Very suspicious | Solid nodule >= 15mm, growing solid/part-solid | Tissue sampling |
| 4X | Very suspicious + | 4A/4B features + additional suspicious morphology | Tissue sampling |
| S | Other | Clinically significant non-pulmonary finding | Per finding type |

**Cross-modal trigger:** Categories 4A, 4B, 4X fire the lung cancer genomics trigger.

**Enum in models.py:**
```python
class LungRADS(str, Enum):
    CAT_0 = "0"
    CAT_1 = "1"
    CAT_2 = "2"
    CAT_3 = "3"
    CAT_4A = "4A"
    CAT_4B = "4B"
    CAT_4X = "4X"
    CAT_S = "S"
```

### C.2 BI-RADS (ACR Breast Imaging)

| Category | Description | Cancer Likelihood | Action |
|---|---|---|---|
| 0 | Incomplete | N/A | Additional imaging |
| 1 | Negative | 0% | Routine screening |
| 2 | Benign | 0% | Routine screening |
| 3 | Probably benign | <2% | Short-interval follow-up |
| 4 | Suspicious | 2-95% | Biopsy |
| 5 | Highly suggestive | >95% | Biopsy |
| 6 | Known malignancy | 100% | Treatment |

### C.3 TI-RADS (ACR Thyroid Imaging)

| Category | Description | FNA Recommendation |
|---|---|---|
| TR1 | Benign | No FNA |
| TR2 | Not suspicious | No FNA |
| TR3 | Mildly suspicious | FNA if >= 2.5 cm; follow if >= 1.5 cm |
| TR4 | Moderately suspicious | FNA if >= 1.5 cm; follow if >= 1.0 cm |
| TR5 | Highly suspicious | FNA if >= 1.0 cm; follow if >= 0.5 cm |

### C.4 LI-RADS (ACR Liver Imaging)

Used for HCC screening in at-risk populations.

| Category | Description | HCC Probability |
|---|---|---|
| LR-1 | Definitely benign | 0% |
| LR-2 | Probably benign | ~5% |
| LR-3 | Intermediate probability | ~30% |
| LR-4 | Probably HCC | ~80% |
| LR-5 | Definitely HCC | ~95% |
| LR-M | Malignant, not HCC specific | Variable |
| LR-TIV | Tumor in vein | Variable |

### C.5 CAD-RADS 2.0 (Coronary Artery Disease)

Used by: `ct_coronary_angiography` workflow

| Category | Stenosis | Description | Management |
|---|---|---|---|
| 0 | 0% | No stenosis | None |
| 1 | 1-24% | Minimal stenosis | Preventive therapy |
| 2 | 25-49% | Mild stenosis | Preventive therapy |
| 3 | 50-69% | Moderate stenosis | Consider functional testing |
| 4A | 70-99% | Severe stenosis | Consider ICA |
| 4B | LM >50% or 3-vessel | Severe stenosis | ICA recommended |
| 5 | 100% | Total occlusion | Consider ICA |
| N | N/A | Non-diagnostic | Repeat or alternative |

**Cross-modal trigger:** Categories 3, 4A, 4B, 5 fire the cardiovascular genomics trigger.

### C.6 PI-RADS v2.1 (Prostate Imaging)

Used by: `mri_prostate_pirads` workflow

| Score | Description | Likelihood of csPC |
|---|---|---|
| 1 | Very low | Very low |
| 2 | Low | Low |
| 3 | Intermediate | Equivocal |
| 4 | High | High |
| 5 | Very high | Very high |

**Cross-modal trigger:** Scores 4 and 5 fire the prostate cancer genomics trigger.

### C.7 ASPECTS (Alberta Stroke Programme Early CT Score)

| Score | Description | Interpretation |
|---|---|---|
| 10 | Normal | No early ischemic changes |
| 7-9 | Mild | Limited early ischemic changes |
| 4-6 | Moderate | Moderate early ischemic changes |
| 0-3 | Severe | Extensive early ischemic changes |

Scoring: Start at 10, subtract 1 point for each affected MCA territory region showing early ischemic change (caudate, lentiform, insula, internal capsule, M1-M6).

### C.8 CT Windowing Presets

Used by: `src/imaging/preview_generator.py`

| Preset | Center (HU) | Width (HU) | Clinical Use |
|---|---|---|---|
| Brain | 40 | 80 | Gray/white matter differentiation |
| Lung | -600 | 1500 | Lung parenchyma, nodules, airways |
| Bone | 400 | 1800 | Fractures, osseous lesions |
| Abdomen | 50 | 400 | Soft tissue organs, masses |
| Soft Tissue | 50 | 350 | Soft tissue detail |

```python
WINDOW_PRESETS = {
    "brain":       {"center": 40,   "width": 80},
    "lung":        {"center": -600, "width": 1500},
    "bone":        {"center": 400,  "width": 1800},
    "abdomen":     {"center": 50,   "width": 400},
    "soft_tissue": {"center": 50,   "width": 350},
}
```

**Windowing formula:**
```
lower = center - width/2
upper = center + width/2
pixel_display = clip(pixel_hu, lower, upper)
pixel_normalized = (pixel_display - lower) / width
```

### C.9 Preview Generation Parameters

| Parameter | Default | Range | Description |
|---|---|---|---|
| `fps` | 8 | 1-30 | Frames per second |
| `max_frames` | 200 | 1-500 | Maximum frame count |
| `format` | `mp4` | mp4, gif | Output format |
| `plane` | `axial` | axial, sagittal, coronal | Slice orientation |

Axis mapping:
```python
_AXIS_MAP = {
    "axial": 2,     # Superior-Inferior
    "sagittal": 0,   # Left-Right
    "coronal": 1,    # Anterior-Posterior
}
```

---

*This guide covers the Imaging Intelligence Agent codebase as of March 2026. For updates, check the repository history and the PROJECT_BIBLE.md in the docs/ directory.*
