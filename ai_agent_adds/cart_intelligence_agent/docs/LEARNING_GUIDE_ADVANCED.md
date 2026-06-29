# Learning Guide -- Advanced

## CAR-T Intelligence Agent: Deep Internals and Extension Guide

**Author:** Adam Jones
**Date:** March 2026
**Codebase Version:** 21,259 lines of Python across 61 files
**Audience:** Experienced developers who want to understand the internals and extend the system

---

## Prerequisites

Before starting this guide, you should have:

1. **Completed the Foundations guide** -- you understand what the CAR-T Intelligence Agent does, how to run it, and how to issue queries through the UI and API.
2. **Python proficiency** -- you are comfortable with Pydantic v2, asyncio, decorators, abstract base classes, and `concurrent.futures`.
3. **Basic ML/NLP concepts** -- you know what embeddings are, what cosine similarity measures, and how retrieval-augmented generation works at a high level.
4. **Vector database basics** -- you understand that Milvus stores high-dimensional vectors and retrieves the nearest neighbors for a query vector.
5. **Development environment** -- you have the repo cloned, dependencies installed, and can run `pytest tests/` successfully (415 tests, all passing).

**Codebase map for reference:**

```
cart_intelligence_agent/
  src/                 # 12,944 lines -- core engine
    rag_engine.py      #   754 lines -- multi-collection RAG
    collections.py     # 1,004 lines -- Milvus schema + manager
    models.py          #   484 lines -- Pydantic models + enums
    knowledge.py       # 2,249 lines -- knowledge graph (3 dictionaries)
    query_expansion.py # 1,592 lines -- 12 expansion maps
    agent.py           #   309 lines -- autonomous agent
    export.py          # 1,487 lines -- Markdown/JSON/PDF export
    metrics.py         #   404 lines -- Prometheus integration
    scheduler.py       #   226 lines -- APScheduler weekly refresh
    ingest/            # ~4,400 lines -- 15 ingest parsers
    utils/             #   pubmed_client.py
  app/                 # 1,162 lines -- Streamlit UI
    cart_ui.py
  api/                 # 1,033 lines -- FastAPI REST API
    main.py
  config/
    settings.py        #   113 lines -- Pydantic BaseSettings
  tests/               # 4,321 lines -- 415 tests
  scripts/             # 1,686 lines -- seed + setup utilities
```

---

## Chapter 1: Deep Dive into the RAG Engine

The RAG engine (`src/rag_engine.py`, 754 lines) is the central nervous system of the agent. Every query -- whether from the Streamlit UI, the FastAPI endpoint, or the autonomous agent -- flows through `CARTRAGEngine`.

### 1.1 The CARTRAGEngine Class

```python
class CARTRAGEngine:
    def __init__(self, collection_manager, embedder, llm_client,
                 knowledge=None, query_expander=None):
        self.collections = collection_manager
        self.embedder = embedder
        self.llm = llm_client
        self.knowledge = knowledge
        self.expander = query_expander
```

Five dependencies are injected at construction time. This is a deliberate design choice: every external service (Milvus, the embedding model, the LLM, the knowledge graph, and the query expansion module) is injected rather than imported directly. This makes the engine fully testable with mocks (see Chapter 9).

### 1.2 The retrieve() Method -- Line by Line

`retrieve()` is the most important method in the entire codebase. Here is the complete execution flow:

**Step 1: Prepare the search text.**

```python
top_k = top_k_per_collection or settings.TOP_K_PER_COLLECTION  # default: 5
start = time.time()

search_text = query.question
if conversation_context:
    search_text = f"{conversation_context}\n\nCurrent question: {query.question}"
```

When conversation memory is active (the UI tracks the last `MAX_CONVERSATION_CONTEXT=3` exchanges), the prior context is prepended to the current question. This changes the embedding to account for conversational continuity -- a follow-up question like "What about its toxicity?" embeds differently when the prior context mentions "Yescarta."

**Step 2: Embed the query.**

```python
query_embedding = self._embed_query(search_text)
```

This calls `_embed_query()`, which prepends the BGE instruction prefix:

```python
def _embed_query(self, text: str):
    prefix = "Represent this sentence for searching relevant passages: "
    return self.embedder.embed_text(prefix + text)
```

The prefix is critical. BGE-small-en-v1.5 is an asymmetric embedding model trained with instruction-prefix pairs. Without the prefix, retrieval quality drops measurably (typically 5-10% lower recall). Documents are embedded without any prefix -- only queries get the prefix.

**Step 3: Build per-collection filter expressions.**

```python
filter_exprs = {}
for coll in collections_to_search:
    parts = []
    cfg = COLLECTION_CONFIG.get(coll, {})
    if query.target_antigen and cfg.get("has_target_antigen"):
        parts.append(f'target_antigen == "{query.target_antigen}"')
    year_field = cfg.get("year_field")
    if year_field:
        if year_min:
            parts.append(f'{year_field} >= {year_min}')
        if year_max:
            parts.append(f'{year_field} <= {year_max}')
    if parts:
        filter_exprs[coll] = " and ".join(parts)
```

This builds Milvus boolean expressions per collection. Not every collection has a `target_antigen` field (e.g., `cart_manufacturing` and `cart_safety` do not -- see `has_target_antigen` in `COLLECTION_CONFIG`). Not every collection has a `year_field` (only `cart_literature` has `year` and `cart_trials` has `start_year`). The filter expressions are passed to Milvus so it applies them during the ANN search, not after.

**Step 4: Parallel search across all collections.**

```python
all_hits = self._search_all_collections(
    query_embedding, collections_to_search, top_k, filter_exprs,
)
```

Inside `_search_all_collections()`, the engine calls `self.collections.search_all()`, which uses `ThreadPoolExecutor` with `max_workers=len(collections)` (11 threads) to search all collections simultaneously. Each thread calls `collection.search()` on Milvus.

**Step 5: Query expansion.**

```python
if self.expander:
    expanded_hits = self._expanded_search(
        query.question, query_embedding, collections_to_search, top_k,
    )
    all_hits.extend(expanded_hits)
```

The expansion system (Chapter 6) generates additional search terms. The key insight in `_expanded_search()` is the bifurcated strategy:

- If an expansion term is a known target antigen (checked against `_KNOWN_ANTIGENS`), it is used as a Milvus field filter: `target_antigen == "CD19"`. This is a precise, metadata-driven search.
- If the expansion term is NOT an antigen (e.g., "cytokine release syndrome"), it is re-embedded and used for a semantic search across all collections. Expansion hits receive a score penalty (0.7x multiplier for semantic, 0.8x for antigen-filtered) to prevent expansion results from overwhelming direct hits.

**Step 6: Merge, deduplicate, and rank.**

```python
hits = self._merge_and_rank(all_hits)
```

```python
def _merge_and_rank(self, hits: List[SearchHit]) -> List[SearchHit]:
    seen = set()
    unique = []
    for hit in hits:
        if hit.id not in seen:
            seen.add(hit.id)
            unique.append(hit)
    unique.sort(key=lambda h: h.score, reverse=True)
    return unique[:30]
```

Dedup is by record ID (`hit.id`). The cap of 30 is a hard limit to prevent prompt bloat -- the LLM context window should not be overwhelmed with low-relevance evidence.

**Step 7: Knowledge graph augmentation.**

```python
knowledge_context = ""
if self.knowledge:
    knowledge_context = self._get_knowledge_context(query.question)
```

`_get_knowledge_context()` is a keyword-matching router that scans the query for mentions of target antigens, toxicities, manufacturing processes, biomarkers, regulatory products, and immunogenicity topics. It calls the appropriate `get_*_context()` function from `knowledge.py` for each match. The knowledge context is injected into the LLM prompt alongside the vector search evidence.

### 1.3 Score Weighting Math

Inside `_search_all_collections()`, each raw Milvus cosine similarity score is adjusted:

```python
weighted_score = raw_score * (1 + weight)
```

Where `weight` comes from `COLLECTION_CONFIG` (sourced from `settings.py`). The default weights are:

| Collection | Weight | Effective Multiplier |
|---|---|---|
| Literature | 0.20 | 1.20 |
| Trials | 0.16 | 1.16 |
| Constructs | 0.10 | 1.10 |
| Assays | 0.09 | 1.09 |
| Safety | 0.08 | 1.08 |
| Biomarkers | 0.08 | 1.08 |
| Manufacturing | 0.07 | 1.07 |
| RealWorld | 0.07 | 1.07 |
| Regulatory | 0.06 | 1.06 |
| Sequences | 0.06 | 1.06 |
| Genomic | 0.04 | 1.04 |

This means a literature hit with a raw score of 0.80 becomes `0.80 * 1.20 = 0.96`, while a genomic hit with the same raw score becomes `0.80 * 1.04 = 0.832`. Literature and trials are intentionally boosted because they contain the most dense, curated information.

### 1.4 Citation Relevance Scoring

Each hit receives a relevance tier:

```python
if raw_score >= settings.CITATION_HIGH_THRESHOLD:     # 0.75
    relevance = "high"
elif raw_score >= settings.CITATION_MEDIUM_THRESHOLD:  # 0.60
    relevance = "medium"
else:
    relevance = "low"
```

These tiers are surfaced in the prompt as `[high relevance]`, `[medium relevance]`, or `[low relevance]` tags, and the prompt instructs the LLM to "prioritize [high relevance] citations."

### 1.5 The CART_SYSTEM_PROMPT

The system prompt covers 12 CAR-T domains and establishes the agent's behavioral rules:

```python
CART_SYSTEM_PROMPT = """You are a CAR-T cell therapy intelligence agent with deep expertise in:

1. **Target Identification** — antigen biology, expression profiling, tumor specificity
2. **CAR Design** — scFv selection, costimulatory domains (CD28 vs 4-1BB), signaling architecture
3. **Vector Engineering** — lentiviral/retroviral production, transduction efficiency, VCN optimization
4. **In Vitro & In Vivo Testing** — cytotoxicity assays, cytokine profiling, mouse models, persistence
5. **Clinical Development** — trial design, response rates, toxicity management (CRS, ICANS)
6. **Manufacturing** — leukapheresis, T-cell expansion, cryopreservation, release testing, CMC
7. **Safety & Pharmacovigilance** — post-market safety signals, REMS, long-term follow-up, FAERS
8. **Biomarkers** — CRS prediction (ferritin, CRP, IL-6), response biomarkers, MRD monitoring
9. **Regulatory Intelligence** — FDA approval pathways, BLA timelines, breakthrough therapy, RMAT, EMA
10. **Molecular Design** — scFv binding affinity, CDR sequences, humanization, nanobodies
11. **Real-World Evidence** — registry outcomes (CIBMTR), community vs academic, special populations
12. **Genomic Evidence** — patient variant data, ClinVar, AlphaMissense pathogenicity
...
```

The prompt explicitly instructs the LLM to cite evidence using clickable markdown links, think cross-functionally, highlight failure modes, and acknowledge uncertainty. This is the personality of the agent.

### 1.6 Prompt Construction

`_build_prompt()` assembles the final prompt from three layers:

1. **Evidence sections** -- grouped by collection, with citation links and relevance tags
2. **Knowledge graph context** -- injected only if `_get_knowledge_context()` returned non-empty text
3. **Instructions** -- the question plus grounding instructions

```python
return (
    f"## Retrieved Evidence\n\n"
    f"{evidence_text}"
    f"{knowledge_text}\n\n"
    f"---\n\n"
    f"## Question\n\n"
    f"{question}\n\n"
    f"Please provide a comprehensive answer grounded in the evidence above. "
    f"Cite sources using the clickable markdown links provided in each evidence item. "
    f"Prioritize [high relevance] citations. "
    f"Consider cross-functional insights across all stages of CAR-T development."
)
```

Each evidence item is formatted as:

```
1. [Literature:PMID 12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/) [high relevance] (score=0.923) CD19 CAR-T therapy achieves...
```

The `_format_citation()` static method generates clickable PubMed URLs for literature hits (when the ID is numeric) and ClinicalTrials.gov URLs for trial hits (when the ID starts with "NCT"). All other collections use the `[Collection:ID]` format.

---

## Chapter 2: Vector Search Internals

### 2.1 How IVF_FLAT Works

The agent uses IVF_FLAT (Inverted File with Flat quantization) as the Milvus index type. Here is what happens at each stage:

**Index building (at collection creation time):**

```python
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}
```

IVF_FLAT partitions the vector space into `nlist=1024` Voronoi cells using k-means clustering. Each vector is assigned to the cell whose centroid is closest. The "FLAT" part means vectors within each cell are stored without compression -- no quantization loss.

**Search (at query time):**

```python
SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}
```

At query time, the system finds the `nprobe=16` closest cluster centroids to the query vector, then performs brute-force comparison against all vectors within those 16 cells. This means the search examines roughly `16/1024 = 1.56%` of the total vectors, providing a massive speedup over brute-force while maintaining high recall.

### 2.2 Why COSINE over L2 or IP

Three common distance metrics for vector search:

| Metric | Formula | Range | Use Case |
|---|---|---|---|
| COSINE | 1 - cos(a,b) | [0, 2] | Normalized similarity, direction matters |
| L2 (Euclidean) | sqrt(sum((a-b)^2)) | [0, inf) | Magnitude-sensitive |
| IP (Inner Product) | sum(a*b) | (-inf, inf) | Pre-normalized vectors |

COSINE was chosen because BGE-small-en-v1.5 embeddings are not pre-normalized. Cosine similarity measures the angle between vectors, ignoring magnitude. This is important for text embeddings where the magnitude can vary based on text length -- a short query should still match a long document chunk. The returned score ranges from 0 (orthogonal, unrelated) to 1 (identical direction, highly similar).

### 2.3 Why nprobe=16

The `nprobe` parameter controls the recall-latency tradeoff:

- **nprobe=1**: Search only the single closest cluster. Very fast (~0.5ms per collection), but recall drops to ~60-70%.
- **nprobe=16**: Search 16 clusters. Moderate latency (~2-5ms per collection), recall ~95-98%.
- **nprobe=1024**: Search all clusters (equivalent to brute force). Recall is 100%, but latency scales linearly with data size.

With 11 collections being searched in parallel, the per-collection latency of 2-5ms means the total search wall-clock time is typically 5-15ms (limited by the slowest collection). The `nprobe=16` value was chosen empirically to keep P99 latency under 50ms while maintaining recall above 95%.

### 2.4 The BGE Embedding Prefix Trick

BGE-small-en-v1.5 is an asymmetric bi-encoder model. It was trained with instruction prefixes for queries but NOT for documents:

```python
# Query embedding (with prefix)
def _embed_query(self, text: str):
    prefix = "Represent this sentence for searching relevant passages: "
    return self.embedder.embed_text(prefix + text)

# Document embedding (no prefix -- done in ingest pipelines)
texts = [record.to_embedding_text() for record in batch]
embeddings = self.embedder.encode(texts)  # No prefix
```

The prefix shifts the query embedding into a region of the vector space that better aligns with relevant document passages. Without the prefix, the model treats the input as a document rather than a query, reducing the quality of retrieval.

### 2.5 How 384 Dimensions Capture Semantics

BGE-small-en-v1.5 maps text into a 384-dimensional vector space. Each dimension encodes a learned semantic feature. The model was trained on 1 billion+ text pairs, learning to place semantically similar texts close together in this space.

For this application:
- "CRS management" and "cytokine release syndrome treatment" have cosine similarity ~0.87
- "CRS management" and "lentiviral transduction efficiency" have cosine similarity ~0.15
- "CD19 CAR-T efficacy in B-ALL" and "tisagenlecleucel response rate in pediatric leukemia" have cosine similarity ~0.82

The 384-dimensional space is sufficient for biomedical text. Larger models (768-dim, 1024-dim) provide marginal improvements but significantly increase storage and latency. At 384 dimensions with float32, each vector consumes 1,536 bytes -- for a collection with 10,000 records, the index is roughly 15MB.

---

## Chapter 3: Adding a New Collection

This chapter walks through a complete worked example: adding a hypothetical `cart_imaging` collection for storing CAR-T cell imaging data (PET-CT, MRI, flow imaging).

### Step 1: Define the Schema in collections.py

Open `src/collections.py` and add a new schema definition alongside the existing ones:

```python
# ── cart_imaging ──────────────────────────────────────────────────────

IMAGING_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Imaging record identifier",
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
        description="Description of imaging finding for embedding",
    ),
    FieldSchema(
        name="modality",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="PET-CT, MRI, flow_imaging, bioluminescence",
    ),
    FieldSchema(
        name="target_antigen",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Target antigen of the CAR-T product imaged",
    ),
    FieldSchema(
        name="timepoint",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="e.g., day 7, day 28, month 3",
    ),
    FieldSchema(
        name="finding",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Key imaging finding",
    ),
    FieldSchema(
        name="product",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="CAR-T product name",
    ),
    FieldSchema(
        name="patient_id",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Anonymized patient identifier",
    ),
]

IMAGING_SCHEMA = CollectionSchema(
    fields=IMAGING_FIELDS,
    description="CAR-T cell imaging and biodistribution data",
)
```

Then register it in `COLLECTION_SCHEMAS`:

```python
COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    # ... existing 11 entries ...
    "cart_imaging": IMAGING_SCHEMA,
}
```

### Step 2: Create the Pydantic Model in models.py

Add a new enum for imaging modality and a new model class:

```python
class ImagingModality(str, Enum):
    PET_CT = "pet_ct"
    MRI = "mri"
    FLOW_IMAGING = "flow_imaging"
    BIOLUMINESCENCE = "bioluminescence"
    SPECT = "spect"


class ImagingRecord(BaseModel):
    """CAR-T cell imaging record -- maps to cart_imaging collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=3000)
    modality: ImagingModality = ImagingModality.PET_CT
    target_antigen: str = Field("", max_length=100)
    timepoint: str = Field("", max_length=50)
    finding: str = Field("", max_length=500)
    product: str = Field("", max_length=200)
    patient_id: str = Field("", max_length=50)

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.modality:
            parts.append(f"Modality: {self.modality.value}")
        if self.finding:
            parts.append(f"Finding: {self.finding}")
        if self.product:
            parts.append(f"Product: {self.product}")
        return " ".join(parts)
```

### Step 3: Register the Model in collections.py

Add the import and the mapping:

```python
from src.models import (
    # ... existing imports ...
    ImagingRecord,
)

COLLECTION_MODELS: Dict[str, type] = {
    # ... existing 11 entries ...
    "cart_imaging": ImagingRecord,
}
```

### Step 4: Add to COLLECTION_CONFIG in rag_engine.py

```python
COLLECTION_CONFIG = {
    # ... existing 11 entries ...
    "cart_imaging": {
        "weight": settings.WEIGHT_IMAGING,
        "label": "Imaging",
        "has_target_antigen": True,
        "year_field": None,
    },
}
```

### Step 5: Add the Weight to settings.py

```python
class CARTSettings(BaseSettings):
    # ... existing weights ...
    WEIGHT_IMAGING: float = 0.05
```

Adjust the other weights so they still sum to approximately 1.0.

### Step 6: Create the Ingest Parser

Create `src/ingest/imaging_parser.py`:

```python
"""Imaging data ingest pipeline for CAR-T Intelligence Agent."""

from typing import Any, List

from loguru import logger
from pydantic import BaseModel

from src.collections import CARTCollectionManager
from src.models import ImagingRecord

from .base import BaseIngestPipeline


class ImagingIngestPipeline(BaseIngestPipeline):
    """Ingest CAR-T imaging data from structured files or APIs."""

    DEFAULT_COLLECTION = "cart_imaging"

    def fetch(self, **kwargs) -> Any:
        """Fetch imaging data from the configured source."""
        # Implement your data fetching logic here
        # Could read from CSV, a PACS API, or a FHIR server
        raise NotImplementedError("Implement fetch() for your data source")

    def parse(self, raw_data: Any) -> List[BaseModel]:
        """Parse raw imaging data into ImagingRecord instances."""
        records = []
        for item in raw_data:
            try:
                record = ImagingRecord(
                    id=item["id"],
                    text_summary=item["description"],
                    modality=item.get("modality", "pet_ct"),
                    target_antigen=item.get("target", ""),
                    timepoint=item.get("timepoint", ""),
                    finding=item.get("finding", ""),
                    product=item.get("product", ""),
                    patient_id=item.get("patient_id", ""),
                )
                records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse imaging record: {e}")
        return records

    def run(self, collection_name=None, batch_size=32, **kwargs) -> int:
        collection_name = collection_name or self.DEFAULT_COLLECTION
        return super().run(collection_name, batch_size, **kwargs)
```

### Step 7: Add Export Format to export.py

In the `_format_evidence_table()` function, add a new branch:

```python
elif collection_name == "Imaging":
    lines.append("| # | ID | Score | Modality | Product | Timepoint | Finding |")
    lines.append("|---|-----|-------|----------|---------|-----------|---------|")
    for i, hit in enumerate(hits[:10], 1):
        m = hit.metadata
        modality = m.get("modality", "")
        product = m.get("product", "")[:20]
        timepoint = m.get("timepoint", "")
        finding = m.get("finding", "")[:40]
        lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {modality} | {product} | {timepoint} | {finding} |")
```

Do the same in `_build_pdf_evidence_table()` for the PDF export.

### Step 8: Add UI Toggle and Badge in cart_ui.py

In the sidebar collection filter section, add `"cart_imaging"` to the list of selectable collections. In the evidence display section, add appropriate badge styling for "Imaging" results.

### Step 9: Add Test Fixtures

In `tests/conftest.py`, update `mock_collection_manager`:

```python
collection_names = [
    # ... existing 10 ...
    "cart_imaging",
]
```

Add a sample imaging hit to `sample_search_hits`:

```python
SearchHit(
    collection="Imaging",
    id="img-pet-001",
    score=0.76,
    text="PET-CT at day 28 shows CAR-T cell trafficking to tumor site.",
    metadata={"modality": "pet_ct", "product": "Yescarta", "timepoint": "day 28"},
),
```

### Step 10: Run Tests

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/cart_intelligence_agent
python -m pytest tests/ -v
```

Verify all 415 existing tests still pass, and write new tests for the imaging model, parser, and export format.

---

## Chapter 4: Building a Custom Ingest Pipeline

### 4.1 The BaseIngestPipeline ABC

All ingest pipelines inherit from `BaseIngestPipeline` in `src/ingest/base.py`. The contract is:

```python
class BaseIngestPipeline(ABC):
    def __init__(self, collection_manager, embedder):
        self.collection_manager = collection_manager
        self.embedder = embedder

    @abstractmethod
    def fetch(self, **kwargs) -> Any:
        """Retrieve raw data from the upstream source."""
        ...

    @abstractmethod
    def parse(self, raw_data: Any) -> List[BaseModel]:
        """Convert raw data into validated Pydantic model instances."""
        ...

    def embed_and_store(self, records, collection_name, batch_size=32) -> int:
        """Embed and insert into Milvus (provided by base class)."""
        ...

    def run(self, collection_name=None, batch_size=32, **fetch_kwargs) -> int:
        """Orchestrate: fetch -> parse -> embed_and_store."""
        ...
```

You implement `fetch()` and `parse()`. The base class provides `embed_and_store()` and `run()`.

### 4.2 The fetch() -> parse() -> embed_and_store() Pattern

The pipeline runs in three strict phases:

1. **fetch()** -- I/O phase. Makes API calls, reads files, or queries databases. Returns raw data in whatever format the source provides (XML, JSON, CSV rows, etc.).

2. **parse()** -- Validation phase. Converts raw data into typed Pydantic model instances. This is where you extract fields, classify records, and enforce data quality. Returns a `List[BaseModel]`.

3. **embed_and_store()** -- Embedding + storage phase. Provided by the base class. For each record, it:
   - Calls `record.to_embedding_text()` to get the text to embed
   - Batches texts (default batch_size=32) and calls `self.embedder.encode(texts)`
   - Converts Pydantic models to dicts via `.model_dump()`
   - Converts Enum values to their string `.value` (Milvus requires plain strings)
   - Truncates UTF-8 strings to respect VARCHAR byte limits
   - Inserts into Milvus via `self.collection_manager.insert_batch()`

### 4.3 Batch Embedding

The base class processes records in batches of 32 by default:

```python
for i in range(0, len(records), batch_size):
    batch = records[i : i + batch_size]
    texts = [record.to_embedding_text() for record in batch]
    embeddings = self.embedder.encode(texts)
    # ... insert batch
```

Batching is important for GPU utilization. BGE-small-en-v1.5 on a DGX Spark can embed ~500 texts/second in batches of 32, but only ~50 texts/second one at a time.

### 4.4 VARCHAR Truncation

Milvus enforces strict byte limits on VARCHAR fields. The base class handles this:

```python
for key, value in record_dict.items():
    if isinstance(value, Enum):
        record_dict[key] = value.value
    elif isinstance(value, str):
        encoded = value.encode("utf-8")
        if len(encoded) > 2990 and key in ("text_chunk", "text_summary"):
            record_dict[key] = encoded[:2990].decode("utf-8", errors="ignore")
        elif len(encoded) > 490 and key in ("title", "name", "known_toxicities"):
            record_dict[key] = encoded[:490].decode("utf-8", errors="ignore")
```

The limits are 2990 bytes (not 3000) and 490 bytes (not 500) to leave a small buffer. UTF-8 characters can be up to 4 bytes, so truncating at byte boundaries with `errors="ignore"` safely drops any incomplete multi-byte characters at the boundary.

### 4.5 Worked Example: Building a FAERS Ingest Pipeline

Here is how the FAERS (FDA Adverse Event Reporting System) pipeline would be structured:

```python
class FAERSIngestPipeline(BaseIngestPipeline):
    DEFAULT_COLLECTION = "cart_safety"

    def fetch(self, product_names=None, start_date=None, **kwargs):
        """Fetch adverse event reports from FAERS API."""
        # Build API query for CAR-T products
        products = product_names or [
            "KYMRIAH", "YESCARTA", "TECARTUS",
            "BREYANZI", "ABECMA", "CARVYKTI"
        ]
        results = []
        for product in products:
            response = requests.get(
                f"https://api.fda.gov/drug/event.json",
                params={
                    "search": f'patient.drug.openfda.brand_name:"{product}"',
                    "limit": 100,
                }
            )
            if response.ok:
                results.extend(response.json().get("results", []))
        return results

    def parse(self, raw_data):
        """Convert FAERS JSON into SafetyRecord instances."""
        records = []
        for event in raw_data:
            record = SafetyRecord(
                id=f"faers-{event['safetyreportid']}",
                text_summary=self._build_summary(event),
                product=self._extract_product(event),
                event_type=self._classify_event(event),
                severity_grade=event.get("serious", "unknown"),
                reporting_source="FAERS",
                year=int(event.get("receiptdate", "0000")[:4]) or 0,
            )
            records.append(record)
        return records

    def run(self, collection_name=None, batch_size=32, **kwargs):
        collection_name = collection_name or self.DEFAULT_COLLECTION
        return super().run(collection_name, batch_size, **kwargs)
```

### 4.6 Error Handling

The base class wraps each batch insertion in a try/except and continues to the next batch on failure:

```python
except Exception as exc:
    logger.error(
        f"Failed batch {i // batch_size + 1} "
        f"({i}-{i + len(batch)}) into '{collection_name}': {exc}"
    )
    continue
```

This means a single corrupt record in a batch will fail that entire batch but not the rest of the pipeline. If you need per-record error isolation, validate records in `parse()` before returning them.

---

## Chapter 5: Extending the Knowledge Graph

The knowledge graph (`src/knowledge.py`, 2,249 lines) contains three curated dictionaries that provide structured, authoritative data to augment vector search results.

### 5.1 The Three Knowledge Dictionaries

| Dictionary | Entries | Key Data |
|---|---|---|
| `CART_TARGETS` | 25 | Protein, UniProt ID, expression, diseases, approved products, resistance, toxicity |
| `CART_TOXICITIES` | 8 | Mechanism, grading, incidence, timing, management, biomarkers, risk factors |
| `CART_MANUFACTURING` | 10 | Process description, parameters, failure modes, release criteria |
| `CART_BIOMARKERS` | 15 | Assay method, clinical cutoff, predictive value, evidence level |
| `CART_REGULATORY` | 6 | Approval dates, pivotal trials, designations, REMS, subsequent approvals |
| `CART_IMMUNOGENICITY` | 6 | ADA incidence, humanization strategies, HLA epitopes, testing paradigms |

### 5.2 How to Add a New Knowledge Dictionary

Suppose you want to add a dictionary for manufacturing equipment/platforms. Follow this pattern:

**Step 1: Define the dictionary.**

```python
# In src/knowledge.py

CART_PLATFORMS: Dict[str, Dict[str, Any]] = {
    "clinimacs_prodigy": {
        "manufacturer": "Miltenyi Biotec",
        "type": "Closed automated system",
        "capabilities": [
            "T-cell activation", "Transduction", "Expansion",
            "Wash", "Formulation", "Sampling"
        ],
        "max_volume": "600 mL",
        "typical_batch_size": "1e7 to 1e10 cells",
        "regulatory_status": "GMP-qualified, CE-marked",
        "products_using": ["Point-of-care CAR-T programs"],
        "advantages": [
            "Fully automated", "Closed system reduces contamination",
            "Single-use tubing set", "Reproducible process"
        ],
        "limitations": [
            "Single batch at a time", "Limited scalability",
            "Operator training required", "High consumable cost"
        ],
    },
    # ... more entries ...
}
```

**Step 2: Write the get_*_context() function.**

```python
def get_platform_context(platform: str) -> str:
    """Return formatted knowledge for a manufacturing platform."""
    key = platform.lower().replace(" ", "_").replace("-", "_")
    data = CART_PLATFORMS.get(key)
    if not data:
        for k in CART_PLATFORMS:
            if key in k.lower():
                data = CART_PLATFORMS[k]
                break
    if not data:
        return ""

    lines = [f"## Platform: {key.replace('_', ' ').title()}"]
    lines.append(f"- **Manufacturer:** {data['manufacturer']}")
    lines.append(f"- **Type:** {data['type']}")
    if data.get("capabilities"):
        lines.append(f"- **Capabilities:** {', '.join(data['capabilities'])}")
    if data.get("advantages"):
        lines.append("- **Advantages:**")
        for a in data["advantages"]:
            lines.append(f"  - {a}")
    return "\n".join(lines)
```

**Step 3: Register in get_all_context_for_query().**

Add a new keyword-matching block:

```python
# Check platforms
platform_keywords = {
    "clinimacs_prodigy": ["PRODIGY", "CLINIMACS"],
    "lonza_cocoon": ["COCOON", "LONZA"],
    # ...
}
for plat_id, keywords in platform_keywords.items():
    if any(kw in query_upper for kw in keywords):
        ctx = get_platform_context(plat_id)
        if ctx:
            sections.append(ctx)
```

**Step 4: Add entity aliases.**

In `ENTITY_ALIASES`:

```python
"PRODIGY": {"type": "platform", "canonical": "clinimacs_prodigy"},
"CLINIMACS": {"type": "platform", "canonical": "clinimacs_prodigy"},
```

**Step 5: Update get_knowledge_stats().**

```python
def get_knowledge_stats() -> Dict[str, int]:
    return {
        # ... existing entries ...
        "manufacturing_platforms": len(CART_PLATFORMS),
    }
```

### 5.3 Keyword Routing

The `get_all_context_for_query()` function uses a simple but effective routing strategy: uppercase string matching against keyword lists. Each knowledge domain has its own block with domain-specific keywords. The function returns early per-domain (using `break` for manufacturing and biomarkers) to avoid injecting too much context from a single domain.

This design is intentionally simple. A future improvement would use the embedding model to classify which knowledge domains are relevant, but keyword matching is fast (sub-millisecond) and reliable for known terms.

### 5.4 Entity Resolution for Comparisons

The `resolve_comparison_entity()` function resolves free-text entity names to structured knowledge entries. It checks five sources in priority order:

1. Target antigens (exact match against `CART_TARGETS` keys)
2. Product/alias table (`ENTITY_ALIASES` -- 54 entries)
3. Toxicities (exact match against `CART_TOXICITIES` keys)
4. Manufacturing processes (fuzzy substring match)
5. Biomarkers (case-insensitive match)

---

## Chapter 6: Query Expansion Engineering

The query expansion system (`src/query_expansion.py`, 1,592 lines) contains 12 expansion maps with hundreds of keyword-to-term mappings.

### 6.1 How Expansion Works

```python
def expand_query(query: str) -> List[str]:
    query_lower = query.lower()
    matched_terms: Set[str] = set()

    for category, mapping in ALL_EXPANSION_MAPS:
        for keyword, terms in mapping.items():
            if keyword in query_lower:
                matched_terms.update(terms)
    return sorted(matched_terms)
```

The function scans the lowercased query for every keyword in every map. When a keyword matches, all of its associated terms are added to a set (automatic deduplication).

### 6.2 The 12 Expansion Maps

```python
ALL_EXPANSION_MAPS: List[tuple] = [
    ("Target Antigen",  TARGET_ANTIGEN_EXPANSION),    # 27 entries
    ("Disease",         DISEASE_EXPANSION),            # 16 entries
    ("Toxicity",        TOXICITY_EXPANSION),           # 12 entries
    ("Manufacturing",   MANUFACTURING_EXPANSION),      # 14 entries
    ("Mechanism",       MECHANISM_EXPANSION),           # 14 entries
    ("Construct",       CONSTRUCT_EXPANSION),           # 18 entries
    ("Safety",          SAFETY_EXPANSION),              #  8 entries
    ("Biomarker",       BIOMARKER_EXPANSION),           # 13 entries
    ("Regulatory",      REGULATORY_EXPANSION),          #  8 entries
    ("Sequence",        SEQUENCE_EXPANSION),            #  8 entries
    ("RealWorld",       REALWORLD_EXPANSION),           # 10 entries
    ("Immunogenicity",  IMMUNOGENICITY_EXPANSION),      # 12 entries
]
```

### 6.3 Writing Effective Expansion Maps

Each map entry is a keyword and a list of semantically related terms:

```python
"crs": [
    "CRS", "cytokine release syndrome", "cytokine storm",
    "tocilizumab", "Actemra", "siltuximab",
    "IL-6", "IL-6 receptor", "sIL-6R",
    "ferritin", "CRP", "C-reactive protein",
    "fever", "hypotension", "hypoxia",
    "Lee grading", "ASTCT grading",
    "grade 3 CRS", "grade 4 CRS",
    "vasopressors", "dexamethasone",
],
```

**Keyword selection strategy:**

1. Use lowercase, short, unambiguous keywords. "crs" is better than "cytokine release syndrome" because it matches more user inputs.
2. Avoid keywords that are common English words. "expansion" as a keyword matches both "T-cell expansion" (manufacturing) and "expansion of clinical trials" (unrelated). If ambiguity is unavoidable, rely on the scoring system to deprioritize irrelevant results.
3. Include product names as keywords in disease maps. Users asking about "multiple myeloma" need to find BCMA-targeted products.

**Term list design:**

1. Include the canonical term first (e.g., "CRS" for the "crs" keyword).
2. Include synonyms and abbreviations ("cytokine release syndrome", "cytokine storm").
3. Include related treatment names ("tocilizumab", "dexamethasone").
4. Include biomarkers ("ferritin", "CRP", "IL-6").
5. Include grading systems ("Lee grading", "ASTCT grading").
6. Keep lists to 10-25 terms. Larger lists increase search breadth but dilute precision.

### 6.4 Avoiding Expansion Explosion

The RAG engine caps expansion to the first 5 expansion terms:

```python
for term in expanded_terms[:5]:
```

Without this cap, a query about "CD19 CAR-T CRS" could expand into 60+ terms, each generating 2+ additional hits across 11 collections -- leading to hundreds of low-relevance results drowning out the primary search.

Expansion hits are also score-penalized:

- Antigen field-filter expansion hits: `score * 0.8`
- Semantic re-embedding expansion hits: `score * 0.7`

This ensures expansion results appear after direct hits in the ranked list.

### 6.5 Category-Aware Expansion

The `expand_query_by_category()` function returns terms grouped by their source map:

```python
>>> expand_query_by_category("transduction efficiency for CD19 CAR")
{
    'Target Antigen': ['CD19', 'B-ALL', 'DLBCL', ...],
    'Manufacturing': ['transduction efficiency', 'lentiviral vector', ...],
}
```

This is useful for weighted category-specific searches where manufacturing terms should receive higher weight when searching `cart_manufacturing`.

### 6.6 Testing Expansion

Use `get_expansion_stats()` to verify your maps are registered:

```python
>>> from src.query_expansion import get_expansion_stats
>>> stats = get_expansion_stats()
>>> for category, data in stats.items():
...     print(f"{category}: {data['keywords']} keywords, {data['total_terms']} total terms")
Target Antigen: 27 keywords, 398 total terms
Disease: 16 keywords, 224 total terms
...
```

---

## Chapter 7: The Comparative Analysis System

### 7.1 How retrieve_comparative() Works

When a user asks "Compare Yescarta vs Kymriah for DLBCL," the system detects this as a comparative query and triggers `retrieve_comparative()`.

**Detection:**

```python
def _is_comparative(self, question: str) -> bool:
    q_upper = question.upper()
    return ("COMPARE" in q_upper or " VS " in q_upper
            or "VERSUS" in q_upper or "COMPARING" in q_upper)
```

**Entity parsing:**

```python
def _parse_comparison_entities(self, question: str):
    # Try "X vs Y" pattern first
    match = re.search(
        r'(.+?)\s+(?:vs\.?|versus)\s+(.+)$',
        q, re.IGNORECASE,
    )
    # Fallback: "compare X and/with Y"
    if not match:
        match = re.search(
            r'(?:compare|comparing)\s+(.+?)\s+(?:and|with)\s+(.+?)(?:\s+(?:for|in)\b.*)?$',
            q, re.IGNORECASE,
        )
```

The parser strips trailing qualifiers ("costimulatory domains", "for DLBCL", "efficacy") and resolves each entity through `resolve_comparison_entity()` in the knowledge graph.

**Dual retrieval:**

```python
evidence_a = self.retrieve(query_a, ...)  # Yescarta-focused search
evidence_b = self.retrieve(query_b, ...)  # Kymriah-focused search
```

Two independent `retrieve()` calls are made, each with the entity's target antigen set as a field filter. This ensures each entity gets its own set of evidence.

**Comparison context:**

```python
if self.knowledge:
    comparison_context = get_comparison_context(entity_a, entity_b)
```

`get_comparison_context()` calls the appropriate `get_*_context()` functions for both entities and joins them with a separator. This provides structured knowledge (approval dates, toxicity profiles, mechanism data) that the LLM can use to build a comparison table.

### 7.2 The Comparative Prompt

The comparative prompt instructs the LLM to produce structured output:

```python
f"1. A **comparison table** in markdown format with key dimensions "
f"as rows and the two entities as columns.\n"
f"2. **Advantages** of each entity (bulleted list).\n"
f"3. **Limitations** of each entity (bulleted list).\n"
f"4. A **clinical context** paragraph explaining when each might "
f"be preferred.\n\n"
```

### 7.3 ComparativeResult Model

The result is a `ComparativeResult` containing two `CrossCollectionResult` instances:

```python
class ComparativeResult(BaseModel):
    query: str
    entity_a: str
    entity_b: str
    evidence_a: CrossCollectionResult
    evidence_b: CrossCollectionResult
    comparison_context: str = ""
    total_search_time_ms: float = 0.0

    @property
    def total_hits(self) -> int:
        return self.evidence_a.hit_count + self.evidence_b.hit_count
```

---

## Chapter 8: Export System Deep Dive

The export system (`src/export.py`, 1,487 lines) generates reports in three formats.

### 8.1 Markdown Export

`export_markdown()` produces a human-readable report with:

- Header: query, timestamp, filters
- Response: the LLM-generated answer verbatim
- Evidence: collection-specific tables using `_format_evidence_table()`
- Knowledge context: verbatim knowledge graph text
- Search metrics: results count, collections searched, latency
- Footer: version string

Each collection has its own table column layout. For example, Literature gets `ID | Score | Source | Title | Year | Target | Journal`, while Safety gets `ID | Score | Product | Event | Severity | Onset | Source`.

### 8.2 JSON Export

`export_json()` produces machine-readable output using Pydantic's `.model_dump()` for proper serialization:

```python
data["evidence"] = evidence.model_dump()
```

This recursively serializes all `SearchHit` objects, including their `metadata` dicts, into a JSON-serializable structure.

### 8.3 PDF Export Architecture

The PDF export uses reportlab's Platypus (Page Layout and Typography Using Scripts) framework. Key components:

**Color palette:**

```python
_NVIDIA_GREEN = colors.HexColor("#76B900")
_DARK_BG = colors.HexColor("#1a1a1a")
_TABLE_ALT = colors.HexColor("#f0f0f0")
_LIGHT_GRAY = colors.HexColor("#666666")
```

**Custom styles:**

```python
def _build_pdf_styles() -> dict:
    return {
        "title": ParagraphStyle("PDFTitle", fontSize=22, textColor=_NVIDIA_GREEN),
        "h2": ParagraphStyle("PDFH2", fontSize=14, textColor=_NVIDIA_GREEN),
        "h3": ParagraphStyle("PDFH3", fontSize=11, textColor="#333333"),
        "body": ParagraphStyle("PDFBody", fontSize=9, leading=13),
        "meta": ParagraphStyle("PDFMeta", fontSize=9, textColor=_LIGHT_GRAY),
        "footer": ParagraphStyle("PDFFooter", fontSize=8, textColor=_LIGHT_GRAY),
    }
```

**Document construction:**

```python
buffer = io.BytesIO()
doc = SimpleDocTemplate(
    buffer, pagesize=letter,
    leftMargin=0.6*inch, rightMargin=0.6*inch,
    topMargin=0.6*inch, bottomMargin=0.6*inch,
)
```

The PDF is built into an in-memory `BytesIO` buffer, then returned as `bytes`. No temporary files are created.

**Evidence tables:**

Each collection gets a styled Table with:
- NVIDIA green header row (`_NVIDIA_GREEN` background, white text)
- Alternating row colors (white and `_TABLE_ALT` gray)
- Grid lines in light gray
- Wrapped text using Paragraph objects in cells (7pt font)

**Markdown-to-flowables conversion:**

The `_md_to_flowables()` function converts the LLM's markdown response into reportlab flowables:

- `## headings` become `Paragraph` with `h2` style
- `### headings` become `Paragraph` with `h3` style
- `**bold**` becomes `<b>bold</b>` XML tags
- `[text](url)` becomes `<a href="url" color="#76B900">text</a>` hyperlinks
- Markdown tables are parsed and rendered as styled reportlab `Table` objects
- Bullet lists (`- item`) are converted to `&#8226; item`
- Block quotes (`> text`) get indented italic styling

### 8.4 Adding a New Export Format

To add an export format (e.g., DOCX), follow this pattern:

1. Create a new public function: `export_docx(query, response_text, evidence, comp_result, filters_applied) -> bytes`
2. Build the document using the same data (query, response, evidence, knowledge context, metrics)
3. Return the content as bytes
4. Add a download button in the Streamlit UI

---

## Chapter 9: Testing Strategies

### 9.1 The Mock-Everything Approach

The test suite runs without Milvus, without the embedding model, and without the LLM. This is achieved through the fixtures in `tests/conftest.py`:

```python
@pytest.fixture
def mock_embedder():
    """384-dim zero vectors."""
    embedder = MagicMock()
    embedder.embed_text.return_value = [0.0] * 384
    return embedder

@pytest.fixture
def mock_llm_client():
    """Always responds 'Mock response'."""
    client = MagicMock()
    client.generate.return_value = "Mock response"
    client.generate_stream.return_value = iter(["Mock ", "response"])
    return client

@pytest.fixture
def mock_collection_manager():
    """search() returns [], search_all() returns {}, stats returns 42."""
    manager = MagicMock()
    manager.search.return_value = []
    manager.search_all.return_value = {name: [] for name in collection_names}
    manager.get_collection_stats.return_value = {name: 42 for name in collection_names}
    return manager
```

The `mock_embedder` returns 384-dimensional zero vectors. This is sufficient for testing the flow -- we are not testing embedding quality, we are testing that the engine correctly assembles, weights, deduplicates, and formats results.

### 9.2 Unit Test Patterns

**Testing the RAG engine:**

```python
def test_retrieve_returns_cross_collection_result(
    mock_embedder, mock_llm_client, mock_collection_manager
):
    engine = CARTRAGEngine(
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
        llm_client=mock_llm_client,
    )
    query = AgentQuery(question="What is CD19?")
    result = engine.retrieve(query)

    assert isinstance(result, CrossCollectionResult)
    assert result.query == "What is CD19?"
    assert result.total_collections_searched == 11
```

**Testing models:**

```python
def test_clinical_trial_id_validation():
    with pytest.raises(ValidationError):
        ClinicalTrial(
            id="INVALID",  # Must match ^NCT\d{8}$
            title="Test",
            text_summary="Test",
        )
```

**Testing knowledge graph:**

```python
def test_get_target_context_cd19():
    ctx = get_target_context("CD19")
    assert "B-Lymphocyte Antigen CD19" in ctx
    assert "ELIANA" in ctx
    assert "Kymriah" in ctx
```

**Testing query expansion:**

```python
def test_crs_expansion_includes_tocilizumab():
    terms = expand_query("What causes CRS?")
    assert "tocilizumab" in terms
    assert "IL-6" in terms
```

**Testing export:**

```python
def test_export_markdown_structure(sample_evidence):
    md = export_markdown(
        query="Test query",
        response_text="Test response",
        evidence=sample_evidence,
    )
    assert "# CAR-T Intelligence Report" in md
    assert "Test query" in md
    assert "Evidence Sources" in md
```

### 9.3 Testing Without Milvus

The entire test suite runs without Milvus because `mock_collection_manager` replaces all Milvus interactions with `MagicMock` return values. When you need to test actual Milvus behavior:

```python
@pytest.fixture(scope="session")
def live_manager():
    """For integration tests only -- requires running Milvus."""
    manager = CARTCollectionManager(host="localhost", port=19530)
    manager.connect()
    yield manager
    manager.disconnect()

@pytest.mark.integration
def test_insert_and_search(live_manager):
    """Requires Milvus. Run with: pytest -m integration"""
    live_manager.create_collection("test_coll", LITERATURE_SCHEMA, drop_existing=True)
    # ... insert records, search, verify results ...
    live_manager.drop_collection("test_coll")
```

Mark integration tests with `@pytest.mark.integration` so they are skipped during normal CI runs.

### 9.4 Property-Based Testing Ideas

For properties that should always hold:

```python
from hypothesis import given, strategies as st

@given(st.text(min_size=1, max_size=500))
def test_expand_query_never_raises(query):
    """expand_query should handle any input without exceptions."""
    result = expand_query(query)
    assert isinstance(result, list)

@given(st.lists(
    st.builds(SearchHit, collection=st.text(), id=st.text(), score=st.floats(0,1), text=st.text()),
    min_size=0, max_size=100,
))
def test_merge_and_rank_caps_at_30(hits):
    """Merged results should never exceed 30."""
    engine = CARTRAGEngine(MagicMock(), MagicMock(), MagicMock())
    result = engine._merge_and_rank(hits)
    assert len(result) <= 30
```

---

## Chapter 10: Performance Optimization

### 10.1 Parallel Search Tuning

The `search_all()` method uses `ThreadPoolExecutor(max_workers=len(collections))` -- 11 threads for 11 collections. This is I/O-bound work (waiting for Milvus network responses), so more threads than CPU cores is appropriate.

If you add more collections and observe thread contention, you can cap the thread pool:

```python
with ThreadPoolExecutor(max_workers=min(len(collections), 20)) as executor:
```

### 10.2 Embedding Caching

For frequently repeated queries (e.g., the same entity lookup in `find_related()`), embedding the same text repeatedly wastes GPU cycles. A simple LRU cache:

```python
from functools import lru_cache

class CachedEmbedder:
    def __init__(self, base_embedder, maxsize=1024):
        self._base = base_embedder
        self._embed = lru_cache(maxsize=maxsize)(self._compute)

    def _compute(self, text: str) -> tuple:
        return tuple(self._base.embed_text(text))

    def embed_text(self, text: str) -> list:
        return list(self._embed(text))
```

Wrap the embedder before injecting it into the engine. The `lru_cache` requires hashable inputs, so we convert the list to a tuple internally.

### 10.3 Milvus Index Parameters

The two critical parameters:

**nlist (index build):** Currently 1024. This creates 1024 cluster centroids. The rule of thumb is `nlist = 4 * sqrt(N)` where N is the number of vectors. For 10,000 vectors: `4 * sqrt(10000) = 400`. For 100,000 vectors: `4 * sqrt(100000) = 1264`. The current value of 1024 is appropriate for collections up to ~65,000 vectors.

**nprobe (search):** Currently 16. Increasing nprobe improves recall at the cost of latency. Benchmark your specific data:

| nprobe | Recall@10 | Latency (ms) |
|---|---|---|
| 4 | ~85% | ~1.5 |
| 16 | ~97% | ~4.0 |
| 64 | ~99% | ~12.0 |
| 256 | ~99.9% | ~40.0 |

### 10.4 Score Threshold Tuning

The default `SCORE_THRESHOLD=0.4` filters out results with cosine similarity below 0.4. This is intentionally permissive -- better to include a marginally relevant result than miss important evidence. If you are seeing too many low-quality results:

```bash
export CART_SCORE_THRESHOLD=0.5
```

Monitor the effect on evidence count per query using the `cart_evidence_count` Prometheus histogram.

### 10.5 Batch Size Optimization

The default `EMBEDDING_BATCH_SIZE=32` balances GPU memory usage and throughput. On the DGX Spark with 128GB unified memory:

| Batch Size | Throughput (texts/sec) | GPU Memory (MB) |
|---|---|---|
| 8 | ~200 | ~100 |
| 32 | ~500 | ~200 |
| 64 | ~600 | ~350 |
| 128 | ~650 | ~600 |

Beyond 64, throughput gains are minimal because BGE-small is a lightweight model. Use 32 as the default; increase to 64 if ingesting large datasets.

---

## Chapter 11: Production Deployment

### 11.1 Docker Multi-Stage Build

Use a multi-stage Dockerfile to minimize image size:

```dockerfile
# Stage 1: Build
FROM python:3.10-slim AS builder
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir --prefix=/install -r requirements.txt

# Stage 2: Runtime
FROM python:3.10-slim
WORKDIR /app
COPY --from=builder /install /usr/local
COPY . .

# Download embedding model at build time
RUN python -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('BAAI/bge-small-en-v1.5')"

EXPOSE 8521 8522
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8522"]
```

The key optimization is downloading the embedding model at build time. This prevents a ~300MB download on every container start.

### 11.2 docker-compose Production Config

```yaml
version: '3.8'
services:
  cart-api:
    build: .
    ports:
      - "8522:8522"
    environment:
      - CART_MILVUS_HOST=milvus
      - CART_MILVUS_PORT=19530
      - CART_LLM_PROVIDER=anthropic
      - ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
      - CART_METRICS_ENABLED=true
      - CART_INGEST_ENABLED=true
    depends_on:
      milvus:
        condition: service_healthy
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8522/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  cart-ui:
    build: .
    command: streamlit run app/cart_ui.py --server.port 8521
    ports:
      - "8521:8521"
    environment:
      - CART_MILVUS_HOST=milvus
      - CART_MILVUS_PORT=19530
      - ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
    depends_on:
      milvus:
        condition: service_healthy

  milvus:
    image: milvusdb/milvus:v2.4-latest
    ports:
      - "19530:19530"
    volumes:
      - milvus_data:/var/lib/milvus
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:9091/healthz"]
      interval: 10s
      timeout: 5s
      retries: 5

volumes:
  milvus_data:
```

### 11.3 Environment Variable Management

All configuration flows through `CARTSettings` (Pydantic BaseSettings with `env_prefix="CART_"`):

| Variable | Default | Description |
|---|---|---|
| `CART_MILVUS_HOST` | localhost | Milvus server hostname |
| `CART_MILVUS_PORT` | 19530 | Milvus server port |
| `CART_LLM_MODEL` | claude-sonnet-4-6 | Anthropic model ID |
| `CART_TOP_K_PER_COLLECTION` | 5 | Max results per collection |
| `CART_SCORE_THRESHOLD` | 0.4 | Min cosine similarity |
| `CART_WEIGHT_LITERATURE` | 0.20 | Literature collection weight |
| `CART_INGEST_SCHEDULE_HOURS` | 168 | Scheduler interval (hours) |
| `CART_INGEST_ENABLED` | false | Enable background ingest |
| `CART_METRICS_ENABLED` | true | Enable Prometheus metrics |
| `CART_CITATION_HIGH_THRESHOLD` | 0.75 | Score for "high relevance" |
| `CART_CITATION_MEDIUM_THRESHOLD` | 0.60 | Score for "medium relevance" |
| `CART_MAX_CONVERSATION_CONTEXT` | 3 | Prior exchanges in memory |

### 11.4 Health Checks

The FastAPI `/health` endpoint returns collection count and total vector count:

```bash
curl http://localhost:8522/health
# {"status":"healthy","collections":11,"total_vectors":1247}
```

Returns HTTP 503 if the engine is not initialized or Milvus is unavailable.

### 11.5 Monitoring with Grafana

The `src/metrics.py` module exposes Prometheus metrics with the `cart_` prefix. Key metrics:

| Metric | Type | Labels | Description |
|---|---|---|---|
| `cart_query_latency_seconds` | Histogram | query_type | Processing time per query |
| `cart_queries_total` | Counter | query_type, status | Total queries by type and status |
| `cart_evidence_count` | Histogram | -- | Evidence items returned per query |
| `cart_collection_hits_total` | Counter | collection | Hits by collection |
| `cart_llm_tokens_total` | Counter | direction | LLM token usage |
| `cart_collection_size` | Gauge | collection | Current records per collection |
| `cart_last_ingest_timestamp` | Gauge | source | Last ingest run timestamp |
| `cart_active_connections` | Gauge | -- | Active API connections |

Create a Grafana dashboard with panels for:

1. Query latency P50/P95/P99 (histogram)
2. Queries per minute by type (counter rate)
3. Evidence count distribution (histogram)
4. Collection sizes over time (gauge)
5. Last ingest freshness (gauge, alert if >168h stale)

---

## Chapter 12: Integration with HCLS AI Factory

### 12.1 The 3-Stage Pipeline

The HCLS AI Factory runs three stages:

1. **Genomics Pipeline** (120-240 min): FASTQ -> VCF via Parabricks (BWA-MEM2, DeepVariant). Produces 11.7M variants.
2. **RAG/Chat Pipeline** (interactive): ClinVar (~2.7M) + AlphaMissense (71M) + Milvus (3.5M vectors) + Claude. Produces the `genomic_evidence` Milvus collection.
3. **Drug Discovery Pipeline** (8-16 min): MolMIM generation + DiffDock docking + RDKit scoring.

The CAR-T Intelligence Agent is a new module that sits alongside Stage 2, reading from the same Milvus instance and adding 10 domain-specific collections.

### 12.2 The Genomic Evidence Bridge

The `genomic_evidence` collection in Milvus is created by the RAG/Chat Pipeline (Stage 2), not by the CAR-T agent. The CAR-T agent treats it as read-only:

```python
# In collections.py
COLLECTION_MODELS: Dict[str, type] = {
    # ... 10 owned collections ...
    "genomic_evidence": None,  # Read-only (no inserts from this agent)
}
```

The schema matches the RAG/Chat Pipeline's output:

```python
GENOMIC_EVIDENCE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=200),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=384),
    FieldSchema(name="chrom", dtype=DataType.VARCHAR, max_length=10),
    FieldSchema(name="pos", dtype=DataType.INT64),
    FieldSchema(name="gene", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="consequence", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="clinical_significance", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="am_pathogenicity", dtype=DataType.FLOAT),
    FieldSchema(name="am_class", dtype=DataType.VARCHAR, max_length=30),
    # ...
]
```

When a user asks "What genomic variants affect CD19 CAR-T response?", the agent searches this collection alongside all 10 CAR-T collections, combining patient-specific genomic data with published literature and clinical trial evidence.

### 12.3 Sharing Milvus with the RAG/Chat Pipeline

Both the CAR-T agent and the RAG/Chat Pipeline connect to the same Milvus instance at port 19530. They use different collection names to avoid conflicts:

- RAG/Chat Pipeline: `clinvar_variants`, `alpha_missense`, `genomic_evidence`
- CAR-T Agent: `cart_literature`, `cart_trials`, `cart_constructs`, etc.

The shared Milvus instance is the integration point. No additional coordination is needed -- Milvus handles concurrent reads from multiple clients.

### 12.4 The VAST AI OS AgentEngine Model

The `CARTIntelligenceAgent` class follows the VAST AI OS AgentEngine pattern:

```
AgentEngine.run()
  |-> Plan          -> search_plan()
  |-> Execute       -> rag_engine.retrieve()
  |-> Reflect       -> evaluate_evidence()
  |-> Report        -> generate_report()
```

This pattern maps directly to VAST AI OS orchestration, where an AgentEngine receives a task, plans its approach, executes subtasks, reflects on results, and produces output.

---

## Chapter 13: Future Architecture

### 13.1 Multi-Agent Systems

The current single-agent architecture can be extended to a multi-agent system where specialized agents collaborate:

- **Literature Agent**: Focused PubMed/PMC search with citation network analysis
- **Clinical Agent**: ClinicalTrials.gov expert with enrollment prediction
- **Manufacturing Agent**: CMC process optimization using historical batch data
- **Safety Agent**: Pharmacovigilance specialist with FAERS signal detection
- **Orchestrator Agent**: Routes queries to the appropriate specialist agent(s)

### 13.2 Graph Databases for Knowledge

The current knowledge graph is a Python dictionary. As the knowledge grows, migrating to a graph database (Neo4j, Amazon Neptune) would enable:

- Relationship traversal: "Find all products targeting the same antigen as Kymriah"
- Path queries: "What connects CRS to CD28 costimulation through IL-6?"
- Graph neural networks for knowledge-graph-augmented retrieval

### 13.3 Fine-Tuned Domain Embeddings

BGE-small-en-v1.5 is a general-purpose embedding model. Fine-tuning on CAR-T literature pairs would improve retrieval quality for domain-specific jargon:

- "scFv" should be closer to "single-chain variable fragment" than to "standard cable variant"
- "4-1BB" should be closer to "CD137 costimulatory" than to "4-1 basketball"
- "Flu/Cy" should be closer to "lymphodepletion conditioning" than to "influenza/cytomegalovirus"

### 13.4 Real-Time Data Streaming

The current scheduler polls PubMed and ClinicalTrials.gov weekly. A streaming architecture would enable:

- Kafka/Pulsar topics for new publications, trial updates, and safety signals
- Change data capture (CDC) from FAERS updates
- Real-time alerts when new evidence matches a user's saved query

### 13.5 VAST AI OS Integration

The VAST AI OS platform provides:

- AgentEngine for stateful agent execution
- ModelEngine for serving BGE-small-en-v1.5 at scale
- DataEngine for Milvus management and data pipelines
- StorageEngine for the DGX Spark's unified memory architecture

Integration would enable the agent to run as a managed service with automatic scaling, health monitoring, and centralized logging.

---

## Appendix A: Complete API Reference

### GET /health

Returns service health with collection and vector counts.

```bash
curl http://localhost:8522/health
```

Response:
```json
{
  "status": "healthy",
  "collections": 11,
  "total_vectors": 1247
}
```

Returns 503 if the engine or Milvus connection is unavailable.

### GET /collections

Returns all collection names and their record counts.

```bash
curl http://localhost:8522/collections
```

Response:
```json
{
  "collections": [
    {"name": "cart_literature", "record_count": 350},
    {"name": "cart_trials", "record_count": 200},
    {"name": "cart_constructs", "record_count": 45},
    {"name": "cart_assays", "record_count": 120},
    {"name": "cart_manufacturing", "record_count": 80},
    {"name": "cart_safety", "record_count": 95},
    {"name": "cart_biomarkers", "record_count": 60},
    {"name": "cart_regulatory", "record_count": 42},
    {"name": "cart_sequences", "record_count": 35},
    {"name": "cart_realworld", "record_count": 70},
    {"name": "genomic_evidence", "record_count": 150}
  ],
  "total": 11
}
```

### POST /query

Full RAG query: retrieve evidence from Milvus, augment with the knowledge graph, and synthesize an LLM response. Requires both the embedding model and LLM client to be available.

```bash
curl -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What are the resistance mechanisms for CD19 CAR-T therapy?",
    "target_antigen": "CD19",
    "year_min": 2020,
    "year_max": 2026
  }'
```

Request body:
| Field | Type | Required | Description |
|---|---|---|---|
| question | string | Yes | Natural-language question |
| target_antigen | string | No | Filter by target antigen (e.g., CD19, BCMA) |
| collections | list[string] | No | Restrict search to specific collections |
| year_min | int | No | Minimum publication year (1990-2030) |
| year_max | int | No | Maximum publication year (1990-2030) |

Response:
```json
{
  "question": "What are the resistance mechanisms for CD19 CAR-T therapy?",
  "answer": "CD19 CAR-T therapy faces several resistance mechanisms...",
  "evidence": [
    {
      "collection": "Literature",
      "id": "35123456",
      "score": 0.912,
      "text": "CD19 loss occurs in 10-20% of patients...",
      "metadata": {"title": "...", "year": 2023}
    }
  ],
  "knowledge_context": "## Target Antigen: CD19\n...",
  "collections_searched": 11,
  "search_time_ms": 34.5
}
```

### POST /search

Evidence-only retrieval (no LLM). Returns evidence snippets without synthesis. Useful for fast lookups or when the client handles its own synthesis.

```bash
curl -X POST http://localhost:8522/search \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Tocilizumab dosing for CRS",
    "collections": ["cart_safety", "cart_biomarkers"]
  }'
```

Response schema is identical to `/query` but without the `answer` field.

### POST /find-related

Cross-collection entity linking. Finds all evidence related to an entity (product name, target antigen, trial ID, etc.) across all 11 collections.

```bash
curl -X POST http://localhost:8522/find-related \
  -H "Content-Type: application/json" \
  -d '{
    "entity": "Yescarta",
    "top_k": 5
  }'
```

Response:
```json
{
  "entity": "Yescarta",
  "results": {
    "cart_literature": [{"collection": "Literature", "id": "...", "score": 0.89, ...}],
    "cart_trials": [{"collection": "Trial", "id": "NCT02348216", "score": 0.95, ...}],
    "cart_constructs": [{"collection": "Construct", "id": "construct-yescarta", ...}],
    "cart_safety": [{"collection": "Safety", "id": "safety-crs-yescarta", ...}],
    "cart_regulatory": [{"collection": "Regulatory", "id": "reg-yescarta-bla", ...}]
  },
  "total_hits": 22
}
```

### GET /knowledge/stats

Returns statistics about the curated knowledge graph.

```bash
curl http://localhost:8522/knowledge/stats
```

Response:
```json
{
  "target_antigens": 34,
  "targets_with_approved_products": 2,
  "toxicity_profiles": 17,
  "manufacturing_processes": 20,
  "biomarkers": 23,
  "regulatory_products": 6
}
```

### GET /metrics

Prometheus-compatible metrics endpoint. Returns counters in Prometheus exposition format.

```bash
curl http://localhost:8522/metrics
```

Response:
```
# HELP cart_api_requests_total Total API requests
# TYPE cart_api_requests_total counter
cart_api_requests_total 147

# HELP cart_api_query_requests_total Total /query requests
# TYPE cart_api_query_requests_total counter
cart_api_query_requests_total 42

# HELP cart_collection_vectors Number of vectors per collection
# TYPE cart_collection_vectors gauge
cart_collection_vectors{collection="cart_literature"} 350
cart_collection_vectors{collection="cart_trials"} 200
...
```

---

## Appendix B: Configuration Reference

All fields from `config/settings.py` (`CARTSettings`):

| Field | Type | Default | Env Variable | Description |
|---|---|---|---|---|
| PROJECT_ROOT | Path | (computed) | -- | Project root directory |
| DATA_DIR | Path | PROJECT_ROOT/data | -- | Data directory |
| CACHE_DIR | Path | DATA_DIR/cache | -- | Cache directory |
| REFERENCE_DIR | Path | DATA_DIR/reference | -- | Reference data directory |
| RAG_PIPELINE_ROOT | Path | /app/rag-chat-pipeline | CART_RAG_PIPELINE_ROOT | Parent RAG pipeline path |
| MILVUS_HOST | str | localhost | CART_MILVUS_HOST | Milvus server hostname |
| MILVUS_PORT | int | 19530 | CART_MILVUS_PORT | Milvus server port |
| COLLECTION_LITERATURE | str | cart_literature | CART_COLLECTION_LITERATURE | Literature collection name |
| COLLECTION_TRIALS | str | cart_trials | CART_COLLECTION_TRIALS | Trials collection name |
| COLLECTION_CONSTRUCTS | str | cart_constructs | CART_COLLECTION_CONSTRUCTS | Constructs collection name |
| COLLECTION_ASSAYS | str | cart_assays | CART_COLLECTION_ASSAYS | Assays collection name |
| COLLECTION_MANUFACTURING | str | cart_manufacturing | CART_COLLECTION_MANUFACTURING | Manufacturing collection name |
| COLLECTION_SAFETY | str | cart_safety | CART_COLLECTION_SAFETY | Safety collection name |
| COLLECTION_BIOMARKERS | str | cart_biomarkers | CART_COLLECTION_BIOMARKERS | Biomarkers collection name |
| COLLECTION_REGULATORY | str | cart_regulatory | CART_COLLECTION_REGULATORY | Regulatory collection name |
| COLLECTION_SEQUENCES | str | cart_sequences | CART_COLLECTION_SEQUENCES | Sequences collection name |
| COLLECTION_REALWORLD | str | cart_realworld | CART_COLLECTION_REALWORLD | Real-world evidence collection name |
| COLLECTION_GENOMIC | str | genomic_evidence | CART_COLLECTION_GENOMIC | Genomic evidence collection (read-only) |
| EMBEDDING_MODEL | str | BAAI/bge-small-en-v1.5 | CART_EMBEDDING_MODEL | HuggingFace model ID |
| EMBEDDING_DIMENSION | int | 384 | CART_EMBEDDING_DIMENSION | Vector dimension |
| EMBEDDING_BATCH_SIZE | int | 32 | CART_EMBEDDING_BATCH_SIZE | Batch size for embedding |
| LLM_PROVIDER | str | anthropic | CART_LLM_PROVIDER | LLM provider |
| LLM_MODEL | str | claude-sonnet-4-6 | CART_LLM_MODEL | LLM model identifier |
| ANTHROPIC_API_KEY | str | None | ANTHROPIC_API_KEY | Anthropic API key |
| TOP_K_PER_COLLECTION | int | 5 | CART_TOP_K_PER_COLLECTION | Max results per collection |
| SCORE_THRESHOLD | float | 0.4 | CART_SCORE_THRESHOLD | Min cosine similarity |
| WEIGHT_LITERATURE | float | 0.20 | CART_WEIGHT_LITERATURE | Literature weight |
| WEIGHT_TRIALS | float | 0.16 | CART_WEIGHT_TRIALS | Trials weight |
| WEIGHT_CONSTRUCTS | float | 0.10 | CART_WEIGHT_CONSTRUCTS | Constructs weight |
| WEIGHT_ASSAYS | float | 0.09 | CART_WEIGHT_ASSAYS | Assays weight |
| WEIGHT_MANUFACTURING | float | 0.07 | CART_WEIGHT_MANUFACTURING | Manufacturing weight |
| WEIGHT_SAFETY | float | 0.08 | CART_WEIGHT_SAFETY | Safety weight |
| WEIGHT_BIOMARKERS | float | 0.08 | CART_WEIGHT_BIOMARKERS | Biomarkers weight |
| WEIGHT_REGULATORY | float | 0.06 | CART_WEIGHT_REGULATORY | Regulatory weight |
| WEIGHT_SEQUENCES | float | 0.06 | CART_WEIGHT_SEQUENCES | Sequences weight |
| WEIGHT_REALWORLD | float | 0.07 | CART_WEIGHT_REALWORLD | Real-world evidence weight |
| WEIGHT_GENOMIC | float | 0.04 | CART_WEIGHT_GENOMIC | Genomic evidence weight |
| NCBI_API_KEY | str | None | CART_NCBI_API_KEY | NCBI API key (optional, increases rate limit) |
| PUBMED_MAX_RESULTS | int | 5000 | CART_PUBMED_MAX_RESULTS | Max PubMed results per query |
| CT_GOV_BASE_URL | str | https://clinicaltrials.gov/api/v2 | CART_CT_GOV_BASE_URL | ClinicalTrials.gov API base URL |
| API_HOST | str | 0.0.0.0 | CART_API_HOST | FastAPI bind address |
| API_PORT | int | 8522 | CART_API_PORT | FastAPI port |
| STREAMLIT_PORT | int | 8521 | CART_STREAMLIT_PORT | Streamlit UI port |
| METRICS_ENABLED | bool | True | CART_METRICS_ENABLED | Enable Prometheus metrics |
| INGEST_SCHEDULE_HOURS | int | 168 | CART_INGEST_SCHEDULE_HOURS | Scheduler interval (hours) |
| INGEST_ENABLED | bool | False | CART_INGEST_ENABLED | Enable background ingest |
| MAX_CONVERSATION_CONTEXT | int | 3 | CART_MAX_CONVERSATION_CONTEXT | Prior exchanges in memory |
| CITATION_HIGH_THRESHOLD | float | 0.75 | CART_CITATION_HIGH_THRESHOLD | Score for high relevance |
| CITATION_MEDIUM_THRESHOLD | float | 0.60 | CART_CITATION_MEDIUM_THRESHOLD | Score for medium relevance |

---

## Appendix C: Metric Reference

### Prometheus Metrics (src/metrics.py)

All metrics use the `cart_` prefix.

**Histograms:**

| Metric | Labels | Buckets | Description |
|---|---|---|---|
| `cart_query_latency_seconds` | query_type | 0.1, 0.5, 1, 2, 5, 10, 30 | Wall-clock processing time per query |
| `cart_evidence_count` | -- | 0, 5, 10, 15, 20, 25, 30 | Evidence items returned per query |

**Counters:**

| Metric | Labels | Description |
|---|---|---|
| `cart_queries_total` | query_type, status | Total queries processed (rag, agent, comparative, entity_link) |
| `cart_collection_hits_total` | collection | Total hits by collection name |
| `cart_llm_tokens_total` | direction | LLM tokens used (prompt, completion) |

**Gauges:**

| Metric | Labels | Description |
|---|---|---|
| `cart_active_connections` | -- | Current active API connections |
| `cart_collection_size` | collection | Current record count per collection |
| `cart_last_ingest_timestamp` | source | Unix timestamp of last ingest run (pubmed, clinical_trials) |

### Helper Functions

```python
from src.metrics import record_query, record_collection_hits, update_collection_sizes, get_metrics_text

# Record a completed query
record_query(query_type="rag", latency=1.23, hit_count=15, status="success")

# Record per-collection hit counts
record_collection_hits({"cart_literature": 5, "cart_trials": 3, "cart_safety": 2})

# Update collection sizes (typically after get_collection_stats())
update_collection_sizes({"cart_literature": 350, "cart_trials": 200, ...})

# Get Prometheus exposition text for /metrics endpoint
text = get_metrics_text()
```

### Graceful Degradation

If `prometheus_client` is not installed, all metric objects become no-op stubs. The application continues to function without metrics:

```python
try:
    from prometheus_client import Counter, Gauge, Histogram, generate_latest
    _PROMETHEUS_AVAILABLE = True
except ImportError:
    _PROMETHEUS_AVAILABLE = False
    # ... no-op stubs ...
```

### API-Level Metrics

The FastAPI `/metrics` endpoint provides basic request counters even without `prometheus_client`:

```
cart_api_requests_total 147
cart_api_query_requests_total 42
cart_api_search_requests_total 89
cart_api_find_related_requests_total 12
cart_api_errors_total 4
cart_collection_vectors{collection="cart_literature"} 350
```

These are simple in-memory Python dict counters incremented on each request. They reset when the process restarts.

---

*Generated by HCLS AI Factory -- CAR-T Intelligence Agent*
*Author: Adam Jones | March 2026*
