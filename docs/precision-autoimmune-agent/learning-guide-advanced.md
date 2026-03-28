# Learning Guide -- Advanced

## Precision Autoimmune Intelligence Agent: Deep Internals and Extension Guide

**Author:** Adam Jones
**Date:** March 2026
**Codebase Version:** 4,283 lines of Python across 10 source files
**Audience:** Experienced developers who want to understand the internals and extend the system

---

## Prerequisites

Before starting this guide, you should have:

1. **Completed the Foundations guide** -- you understand what the Precision Autoimmune Intelligence Agent does, how to run it, and how to issue queries through the UI and API.
2. **Python proficiency** -- you are comfortable with Pydantic v2, asyncio, decorators, abstract base classes, and `concurrent.futures`.
3. **Basic ML/NLP concepts** -- you know what embeddings are, what cosine similarity measures, and how retrieval-augmented generation works at a high level.
4. **Vector database basics** -- you understand that Milvus stores high-dimensional vectors and retrieves the nearest neighbors for a query vector.
5. **Development environment** -- you have the repo cloned, dependencies installed, and can run `pytest tests/` successfully (455 tests, all passing).

**Codebase map for reference:**

```
precision_autoimmune_agent/
  src/                       # 4,283 lines -- core engine
    agent.py                 #   437 lines -- AutoimmuneAgent orchestrator
    rag_engine.py            #   597 lines -- AutoimmuneRAGEngine, SearchHit, CrossCollectionResult
    collections.py           #   562 lines -- AutoimmuneCollectionManager (14 collections)
    models.py                #   238 lines -- 17 Pydantic models + 3 enums
    knowledge.py             #   855 lines -- 5 knowledge domains
    diagnostic_engine.py     #   519 lines -- DiagnosticEngine + CLASSIFICATION_CRITERIA
    document_processor.py    #   435 lines -- DocumentProcessor (PDF ingestion)
    export.py                #   389 lines -- Markdown/PDF/FHIR R4 export
    timeline_builder.py      #   251 lines -- TimelineBuilder
  app/
    autoimmune_ui.py         # Streamlit UI (port 8531)
  api/
    main.py                  # FastAPI REST API (port 8532)
  config/
    settings.py              # AutoimmuneSettings (Pydantic BaseSettings, AUTO_ prefix)
  tests/                     # 455 tests across 8 files
```

---

## Chapter 1: Deep Dive into the RAG Engine

The RAG engine (`src/rag_engine.py`, 597 lines) is the central nervous system of the agent. Every query -- whether from the Streamlit UI, the FastAPI endpoint, or the diagnostic engine -- flows through `AutoimmuneRAGEngine`.

### 1.1 The AutoimmuneRAGEngine Class

```python
class AutoimmuneRAGEngine:
    def __init__(
        self,
        collection_manager,
        embedder,
        llm_client,
        settings=None,
        knowledge=None,
    ):
        self.cm = collection_manager
        self.embedder = embedder
        self.llm = llm_client
        self.knowledge = knowledge

        if settings is None:
            from config.settings import settings as _settings
            settings = _settings
        self.settings = settings
        self._conversation_lock = threading.Lock()
        self._conversation_history: deque = deque(maxlen=self.settings.CONVERSATION_MEMORY_SIZE)
        self._embed_cache: Dict[str, List[float]] = {}
        self._embed_cache_max = 256
```

Five dependencies are injected at construction time. This is a deliberate design choice: every external service (Milvus, the embedding model, the LLM, the knowledge base, and the settings) is injected rather than imported directly. This makes the engine fully testable with mocks (see Chapter 10).

The settings fallback is noteworthy: if no settings object is provided, the engine lazily imports the global singleton. This means tests can inject a custom settings object while production code gets the environment-configured default.

Two additional internal structures are initialized:

- **`_conversation_history`**: A thread-safe `deque` with a max length of `CONVERSATION_MEMORY_SIZE` (default 3). This stores the most recent Q&A pairs for conversational continuity.
- **`_embed_cache`**: A simple dictionary that functions as a 256-entry LRU cache for query embeddings. This avoids re-encoding identical queries within a session.

### 1.2 The retrieve() Method -- Line by Line

`retrieve()` is the most important method in the entire codebase. Here is the complete execution flow:

**Step 1: Guard on embedder availability.**

```python
if not self.can_search:
    logger.warning("RAG retrieve called but embedder not available")
    return CrossCollectionResult(query=query, search_time_ms=0)
```

The `can_search` property checks `self.embedder is not None`. If the embedding model failed to load at startup (common in memory-constrained environments), the engine returns an empty result rather than crashing. This graceful degradation pattern runs throughout the codebase.

**Step 2: Set top_k from settings if not provided.**

```python
if top_k_per_collection is None:
    top_k_per_collection = self.settings.TOP_K_PER_COLLECTION  # default: 5
```

This means each of the 14 collections can return up to 5 hits, yielding a theoretical maximum of 70 raw hits before deduplication.

**Step 3: Embed the query.**

```python
try:
    query_embedding = self._embed(query)
except Exception as exc:
    logger.error(f"Embedding failed: {exc}")
    return CrossCollectionResult(query=query, search_time_ms=0)
```

This calls `_embed()`, which prepends the BGE instruction prefix:

```python
def _embed(self, text: str) -> List[float]:
    if self.embedder is None:
        raise RuntimeError("Embedder not initialized")
    cache_key = text[:512]
    if cache_key in self._embed_cache:
        return self._embed_cache[cache_key]
    instruction = self.settings.BGE_INSTRUCTION
    vec = self.embedder.encode(instruction + text).tolist()
    # Evict oldest if cache full
    if len(self._embed_cache) >= self._embed_cache_max:
        oldest_key = next(iter(self._embed_cache))
        del self._embed_cache[oldest_key]
    self._embed_cache[cache_key] = vec
    return vec
```

The instruction prefix is critical. BGE-small-en-v1.5 is an asymmetric embedding model trained with instruction-prefix pairs. The prefix is `"Represent this sentence for searching relevant passages: "`. Without it, retrieval quality drops measurably (typically 5-10% lower recall). Documents are embedded without any prefix -- only queries get it.

The cache key truncates to 512 characters. This means two queries that differ only after character 512 will collide. In practice, autoimmune clinical queries rarely exceed this length, so the tradeoff favors cache hit rate.

**Step 4: Detect disease areas and build filter expressions.**

```python
disease_areas = self._detect_disease_areas(query)

filter_exprs = {}
if patient_id:
    safe_pid = _sanitize_filter_value(patient_id)
    if safe_pid:
        for coll_name in (self.settings.COLL_CLINICAL_DOCUMENTS,
                          self.settings.COLL_PATIENT_LABS,
                          self.settings.COLL_PATIENT_TIMELINES):
            filter_exprs[coll_name] = f'patient_id == "{safe_pid}"'
```

Disease area detection uses keyword matching against `DISEASE_KEYWORDS` -- a dictionary of 12 disease areas with 4-10 keywords each (rheumatoid arthritis, systemic lupus, multiple sclerosis, ankylosing spondylitis, Sjogren's, systemic sclerosis, IBD, psoriasis, myasthenia gravis, celiac, thyroid autoimmune, and POTS/EDS/MCAS).

Filter expressions are only applied to patient-scoped collections (clinical documents, labs, timelines). The `_sanitize_filter_value()` function rejects values containing characters outside `[A-Za-z0-9 _\-\.]` and limits length to 64 characters. This prevents Milvus expression injection.

**Step 5: Parallel search across all collections.**

```python
raw_results = self.cm.search_all(
    query_embedding=query_embedding,
    top_k_per_collection=top_k_per_collection,
    collections=collections_filter,
    filter_exprs=filter_exprs,
    score_threshold=self.settings.SCORE_THRESHOLD,  # default: 0.40
)
```

Inside `search_all()`, the collection manager uses `ThreadPoolExecutor` with `max_workers=6` to search collections in parallel. Each thread calls `collection.search()` on Milvus. Results below `SCORE_THRESHOLD` (0.40) are filtered out server-side.

**Step 6: Convert, deduplicate, weight, and rank.**

```python
for coll_name, coll_hits in raw_results.items():
    weight = config.get(coll_name, {}).get("weight", 0.05)
    for h in coll_hits:
        hit_id = h["id"]
        if hit_id in seen_ids:
            continue
        seen_ids.add(hit_id)

        text = h.get("text_chunk", h.get("text_summary", ""))
        text_hash = hashlib.md5(text[:300].encode()).hexdigest()
        if text_hash in seen_texts:
            continue
        seen_texts.add(text_hash)

        weighted_score = min(h["score"] * (1 + weight), 1.0)
        relevance = self._score_relevance(h["score"])
```

Deduplication happens at two levels: by ID and by content hash (MD5 of the first 300 characters). Content-hash deduplication catches cases where the same text appears in different collections (e.g., a lab result mentioned in both `autoimmune_patient_labs` and `autoimmune_clinical_documents`).

**Step 7: Sort and cap.**

```python
hits.sort(key=lambda x: x.score, reverse=True)
hits = hits[: self.settings.MAX_EVIDENCE_ITEMS]  # default: 30
```

After deduplication, results are sorted by raw score (not weighted score -- the weight only affects the final score reported in the `SearchHit`). The top 30 hits are retained to keep the LLM context window manageable.

### 1.3 Score Weighting Math

The weighted score formula is:

```
weighted_score = min(cosine_similarity * (1 + collection_weight), 1.0)
```

This is a multiplicative boost, not an additive one. A hit with cosine similarity 0.85 from `autoimmune_clinical_documents` (weight 0.18) becomes:

```
0.85 * (1 + 0.18) = 0.85 * 1.18 = 1.003 → capped at 1.0
```

The cap at 1.0 prevents inflated scores. The 14 collection weights sum to ~1.0 and are configured in `settings.py`:

| Collection | Weight | Label |
|------------|--------|-------|
| autoimmune_clinical_documents | 0.18 | Clinical Document |
| autoimmune_patient_labs | 0.14 | Lab Result |
| autoimmune_autoantibody_panels | 0.12 | Autoantibody |
| autoimmune_hla_associations | 0.08 | HLA Association |
| autoimmune_disease_criteria | 0.08 | Classification Criteria |
| autoimmune_disease_activity | 0.07 | Disease Activity |
| autoimmune_flare_patterns | 0.06 | Flare Pattern |
| autoimmune_biologic_therapies | 0.06 | Biologic Therapy |
| autoimmune_clinical_trials | 0.05 | Clinical Trial |
| autoimmune_literature | 0.05 | Literature |
| autoimmune_pgx_rules | 0.04 | PGx Rule |
| autoimmune_patient_timelines | 0.03 | Timeline |
| autoimmune_cross_disease | 0.02 | Cross-Disease |
| genomic_evidence | 0.02 | Genomic Evidence |

The weighting reflects clinical priority: patient-specific data (documents, labs, antibodies) is weighted highest because it is most relevant to individual patient queries. Reference data (PGx rules, cross-disease patterns, genomic evidence) is weighted lowest because it provides general context.

### 1.4 Citation Relevance Scoring

The `_score_relevance()` method maps raw cosine similarity to three tiers:

```python
def _score_relevance(self, score: float) -> str:
    if score >= self.settings.CITATION_HIGH:     # >= 0.80
        return "high"
    elif score >= self.settings.CITATION_MEDIUM:  # >= 0.60
        return "medium"
    return "low"
```

The thresholds are configurable via `AUTO_CITATION_HIGH` and `AUTO_CITATION_MEDIUM` environment variables. The relevance tag is surfaced in evidence blocks sent to the LLM and in the UI as colored badges (green/yellow/red).

These thresholds were calibrated empirically. With BGE-small-en-v1.5 on autoimmune clinical text:
- **>= 0.80**: Near-exact semantic match. The evidence directly answers the query.
- **>= 0.60**: Topically relevant. The evidence is about the right disease/test but may not directly answer.
- **< 0.60**: Weak relevance. Included for coverage but should be treated with caution.

### 1.5 The System Prompt

The system prompt (`SYSTEM_PROMPT` in `rag_engine.py`) is a 1,500-character instruction that establishes the agent's clinical identity and output format. Key elements:

1. **Domain expertise declaration**: Autoantibody interpretation, HLA-disease associations (with specific examples like HLA-B*27:05 OR=87.4), disease activity scoring, flare prediction, biologic therapy selection, diagnostic odyssey analysis, ACR/EULAR classification criteria, overlap syndromes, and pharmacogenomics.

2. **Citation format specification**: Six citation types are defined:
   - `[AutoAb:name]` for autoantibodies
   - `[HLA:allele]` for HLA associations
   - `[Activity:score_name]` for activity scores
   - `[Therapy:drug]` for biologics
   - `[Literature:PMID](url)` for published literature
   - `[Trial:NCT_ID](url)` for clinical trials

3. **Clinical safety guardrail**: The prompt ends with `"This is a clinical decision-support tool. All recommendations should be reviewed by a qualified healthcare provider before clinical action."`

### 1.6 Knowledge Augmentation Pipeline

The `_build_knowledge_context()` method enriches the LLM prompt with structured knowledge from four domains:

```python
def _build_knowledge_context(self, query: str, disease_areas: List[str]) -> str:
```

It scans the query text for mentions of:

1. **HLA alleles** -- matches against `HLA_DISEASE_ASSOCIATIONS` (22 alleles). Adds disease name, odds ratio, and PMID.
2. **Autoantibodies** -- matches against `AUTOANTIBODY_DISEASE_MAP` (24 antibodies). Adds sensitivity, specificity, and clinical notes.
3. **Biologic therapies** -- matches drug names and drug classes against `BIOLOGIC_THERAPIES` (22 drugs). Adds mechanism, indications, and PGx considerations.
4. **Flare patterns** -- matches disease names against `FLARE_BIOMARKER_PATTERNS` (13 diseases). Adds early warning biomarkers and protective signals.

Results are capped at 25 knowledge items (`parts[:25]`). This cap prevents knowledge context from overwhelming the retrieved evidence in the LLM prompt.

### 1.7 Conversation Memory (3 turns)

The engine maintains a sliding window of the last 3 conversation exchanges:

```python
self._conversation_history: deque = deque(maxlen=self.settings.CONVERSATION_MEMORY_SIZE)
```

When building messages for the LLM, prior exchanges are prepended with truncated content:

```python
for entry in history:
    messages.append({"role": "user", "content": entry["question"][:200]})
    messages.append({"role": "assistant", "content": entry["answer"][:800]})
```

Questions are truncated to 200 characters and answers to 800 characters. This keeps conversation context concise while providing enough continuity for follow-up questions like "What about the PGx considerations?" to resolve their referent.

Thread safety is ensured with `self._conversation_lock` (a `threading.Lock`). All reads and writes to `_conversation_history` happen inside `with self._conversation_lock:` blocks.

---

## Chapter 2: Vector Search Internals

### 2.1 How IVF_FLAT Works

The agent uses IVF_FLAT (Inverted File with Flat storage) as its Milvus index type. Defined in `collections.py`:

```python
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}
```

IVF_FLAT partitions the vector space into `nlist=1024` Voronoi cells (clusters). At index time, each vector is assigned to its nearest cluster centroid. At search time, only `nprobe` clusters are searched instead of all 1024.

This is the right choice for the autoimmune agent because:
- **Exact within-cluster results**: Unlike IVF_PQ (product quantization) or HNSW, IVF_FLAT stores the original vectors, not compressed approximations. This means no recall loss within searched clusters.
- **Collection sizes**: Most autoimmune collections have 1K-50K vectors. At this scale, IVF_FLAT provides excellent recall with low latency.
- **Predictable memory**: Each vector is 384 floats x 4 bytes = 1,536 bytes. A 50K-vector collection uses ~73MB of index memory.

### 2.2 Why COSINE over L2 or IP

The metric type is COSINE (cosine similarity), not L2 (Euclidean) or IP (inner product):

```python
SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}
```

Cosine similarity measures the angle between two vectors, ignoring magnitude. This is the correct choice for BGE embeddings because:

- **BGE embeddings are not normalized**: While some embedding models output unit vectors (where cosine = IP), BGE-small-en-v1.5 does not guarantee this. Using IP directly would give incorrect rankings.
- **Length-invariant**: A 50-word chunk and a 500-word chunk about the same topic should score similarly. COSINE achieves this; L2 would penalize magnitude differences.
- **Score interpretation**: COSINE scores range from -1 to 1, where 1.0 = identical direction. This makes the threshold at 0.40 and the HIGH/MEDIUM/LOW buckets interpretable.

### 2.3 Why nprobe=16

The `nprobe=16` parameter means that at search time, Milvus scans 16 of the 1024 clusters (~1.56% of the index). This is a latency/recall tradeoff:

| nprobe | Approx recall | Latency impact |
|--------|---------------|----------------|
| 1 | ~65% | Fastest |
| 8 | ~90% | Low |
| **16** | **~95%** | **Moderate** |
| 32 | ~98% | Higher |
| 1024 | 100% (brute force) | Highest |

The value 16 was chosen because it provides > 95% recall while keeping per-collection search latency under 5ms on DGX Spark. Since the agent searches 14 collections in parallel, total search latency is dominated by the slowest collection, not the sum.

### 2.4 The BGE Embedding Prefix Trick

BGE-small-en-v1.5 uses an asymmetric retrieval paradigm with instruction prefixes. The agent uses it as follows:

**Query embedding** (in `_embed()`):
```python
instruction = self.settings.BGE_INSTRUCTION
# "Represent this sentence for searching relevant passages: "
vec = self.embedder.encode(instruction + text).tolist()
```

**Document embedding** (in `DocumentProcessor.embed_records()`):
```python
texts = [r["text_chunk"] for r in records]
embeddings = self.embedder.encode(texts, batch_size=32, show_progress_bar=False)
```

Documents are embedded without any prefix. This asymmetry is by design: the instruction tells the model "I'm searching for something" versus "I am a document to be found." Omitting the prefix from documents during ingestion is correct behavior, not a bug.

### 2.5 How 384 Dimensions Capture Semantics

BGE-small-en-v1.5 produces 384-dimensional vectors. The dimension constant is defined in `collections.py`:

```python
_DIM = 384  # BGE-small-en-v1.5
```

Each dimension encodes a learned semantic feature. Nearby vectors in this 384-dimensional space represent semantically similar concepts. The model was trained on a massive corpus of text pairs, learning that "rheumatoid arthritis joint swelling" should be close to "RA synovitis with morning stiffness."

384 dimensions is a deliberate tradeoff:
- **BGE-large**: 1024 dimensions, 335M parameters. Higher quality but 2.7x more memory per vector.
- **BGE-small**: 384 dimensions, 33M parameters. Good quality with much lower memory footprint.
- For the autoimmune agent's domain-specific vocabulary, the quality difference is marginal because retrieval is supplemented by knowledge augmentation.

### 2.6 Parallel Collection Search (ThreadPoolExecutor)

The `search_all()` method in `AutoimmuneCollectionManager` searches all collections simultaneously:

```python
def search_all(
    self,
    query_embedding: List[float],
    top_k_per_collection: int = 5,
    collections: Optional[List[str]] = None,
    filter_exprs: Optional[Dict[str, str]] = None,
    score_threshold: float = 0.0,
    max_workers: int = 6,
) -> Dict[str, List[Dict[str, Any]]]:
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                self.search, coll_name, query_embedding,
                top_k_per_collection, filter_exprs.get(coll_name),
            ): coll_name
            for coll_name in collections
        }
        for future in as_completed(futures):
            coll_name = futures[future]
            hits = future.result()
            if score_threshold > 0:
                hits = [h for h in hits if h["score"] >= score_threshold]
            if hits:
                results[coll_name] = hits
```

Key design decisions:

- **`max_workers=6`**: Limits concurrent Milvus connections. The DGX Spark has ample CPU threads, but Milvus server-side query execution benefits from limiting client-side concurrency to avoid lock contention.
- **`as_completed()`**: Results are processed as they arrive, not in submission order. The fastest collections contribute hits immediately.
- **Score threshold filtering**: Applied client-side after each collection search returns. This is a post-filter, not a Milvus-level filter.

---

## Chapter 3: Adding a New Collection

Adding a new Milvus collection requires changes across 10 files. Here is the step-by-step process:

### Step 1: Define Schema in collections.py

Add the collection schema using the helper functions:

```python
# 15. Drug interactions
_register(
    "autoimmune_drug_interactions",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("drug_a", 128),
        _varchar("drug_b", 128),
        _varchar("interaction_type", 64),      # synergistic, antagonistic, contraindicated
        _varchar("severity", 32),               # mild, moderate, severe
        _varchar("mechanism", 1024),
        _varchar("clinical_recommendation", 2000),
        _varchar("evidence_level", 64),
        _varchar("pmid", 256),
    ],
    "Drug-drug interaction database for autoimmune therapies",
)
```

Every collection must have `_pk()` (primary key) and `_embedding()` (384-dim float vector) as its first two fields. The `text_chunk` field is required for the RAG engine to extract display text.

### Step 2: Create Pydantic Model in models.py

```python
class DrugInteractionRecord(BaseModel):
    """Record for autoimmune_drug_interactions collection."""
    id: str
    text_chunk: str
    drug_a: str = ""
    drug_b: str = ""
    interaction_type: str = ""
    severity: str = ""
    mechanism: str = ""
    clinical_recommendation: str = ""
    evidence_level: str = ""
    pmid: str = ""

    def to_embedding_text(self) -> str:
        return (
            f"Drug interaction: {self.drug_a} + {self.drug_b} ({self.interaction_type}). "
            f"Severity: {self.severity}. {self.mechanism}. {self.text_chunk}"
        )
```

The `to_embedding_text()` method generates the text that will be embedded. It should front-load the most discriminative content (drug names, interaction type) before appending the raw text chunk.

### Step 3: Register in collections.py

The `_register()` call in Step 1 handles this. The schema is added to `COLLECTION_SCHEMAS`, a module-level dictionary:

```python
COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {}

def _register(name: str, fields: list, description: str = ""):
    schema = CollectionSchema(fields=fields, description=description)
    COLLECTION_SCHEMAS[name] = schema
```

### Step 4: Add Weight to settings.py

Add a weight variable with the `AUTO_` prefix:

```python
class AutoimmuneSettings(BaseSettings):
    # ... existing weights ...
    WEIGHT_DRUG_INTERACTIONS: float = 0.03
```

Add the corresponding collection name constant:

```python
COLL_DRUG_INTERACTIONS: str = "autoimmune_drug_interactions"
```

Update the `_validate_weights` model validator to include the new weight in the sum check. Update the `collection_config` property:

```python
self.COLL_DRUG_INTERACTIONS: {
    "weight": self.WEIGHT_DRUG_INTERACTIONS,
    "label": "Drug Interaction",
    "name": self.COLL_DRUG_INTERACTIONS,
},
```

### Step 5: Add to COLLECTION_CONFIG in rag_engine.py

The RAG engine reads collection config from `self.settings.collection_config`, so after updating `settings.py`, the engine automatically picks up the new collection. No changes needed in `rag_engine.py` itself.

### Step 6: Create Ingest Parser

Create a function that transforms raw data into collection-ready records:

```python
def parse_drug_interactions(source_file: Path) -> List[Dict[str, Any]]:
    """Parse drug interaction data from a TSV/JSON file."""
    records = []
    data = json.loads(source_file.read_text())
    for entry in data:
        record_id = hashlib.md5(
            f"{entry['drug_a']}_{entry['drug_b']}".encode()
        ).hexdigest()[:16]
        records.append({
            "id": f"ddi_{record_id}",
            "text_chunk": entry.get("description", "")[:3000],
            "drug_a": entry["drug_a"],
            "drug_b": entry["drug_b"],
            "interaction_type": entry.get("type", ""),
            "severity": entry.get("severity", ""),
            "mechanism": entry.get("mechanism", ""),
            "clinical_recommendation": entry.get("recommendation", ""),
            "evidence_level": entry.get("evidence", ""),
            "pmid": entry.get("pmid", ""),
        })
    return records
```

### Step 7: Add Export Format to export.py

In `AutoimmuneExporter._format_analysis_md()`, add a section for the new data type if it appears in analysis results. In `to_fhir_r4()`, map the interaction to an appropriate FHIR resource (e.g., `DetectedIssue`).

### Step 8: Add UI Toggle in autoimmune_ui.py

In the Streamlit sidebar, add a checkbox for the new collection:

```python
include_ddi = st.sidebar.checkbox("Drug Interactions", value=True)
if not include_ddi:
    excluded_collections.append("autoimmune_drug_interactions")
```

### Step 9: Add Test Fixtures

Add to the test suite:

```python
@pytest.fixture
def sample_drug_interaction():
    return {
        "id": "ddi_test_001",
        "text_chunk": "Methotrexate and adalimumab combination therapy...",
        "drug_a": "Methotrexate",
        "drug_b": "Adalimumab",
        "interaction_type": "synergistic",
        "severity": "mild",
        "embedding": [0.0] * 384,
    }
```

### Step 10: Run Tests

```bash
pytest tests/ -v --tb=short
# All 455+ tests should pass, including the new ones
```

---

## Chapter 4: Building a Custom Ingest Pipeline

### 4.1 The DocumentProcessor Class

The `DocumentProcessor` (`src/document_processor.py`, 435 lines) implements the full ingestion pipeline from PDF files to Milvus vectors:

```python
class DocumentProcessor:
    def __init__(
        self,
        collection_manager=None,
        embedder=None,
        max_chunk_size: int = 2500,
        chunk_overlap: int = 200,
    ):
```

Four parameters control the processor:
- **`collection_manager`**: The `AutoimmuneCollectionManager` instance for inserting records.
- **`embedder`**: A `SentenceTransformer` instance for generating embeddings.
- **`max_chunk_size`**: Maximum characters per chunk (default 2500). This maps to roughly 500-600 tokens for BGE.
- **`chunk_overlap`**: Characters of overlap between consecutive chunks (default 200). This ensures entities near chunk boundaries appear in at least one chunk.

### 4.2 PDF -> Text -> Chunks -> Entities -> Embeddings -> Milvus

The full pipeline in `process_pdf()`:

```python
def process_pdf(self, pdf_path: Path, patient_id=None) -> List[Dict[str, Any]]:
    pages = self.extract_pages_from_pdf(pdf_path)       # PyPDF2 extraction
    full_text = "\n\n".join(text for _, text in pages)
    doc_type = self.classify_document_type(full_text)    # 7 doc types
    specialty = self.detect_specialty(full_text)          # 11 specialties
    visit_date = self.extract_date(full_text)            # 4 regex patterns
    provider = self.extract_provider(full_text)          # Name extraction

    records = []
    for page_num, page_text in pages:
        chunks = self.chunk_text(page_text)              # Sentence-boundary chunking
        for chunk_idx, chunk in enumerate(chunks):
            record_id = hashlib.md5(
                f"{pdf_path.name}:{page_num}:{chunk_idx}".encode()
            ).hexdigest()[:16]
            records.append({...})
    return records
```

The pipeline is deterministic: given the same PDF, it always produces the same record IDs (via MD5 of filename:page:chunk_index). This makes re-ingestion idempotent -- inserting the same PDF twice creates records with the same IDs that Milvus can detect as duplicates.

### 4.3 Entity Extraction

The document processor extracts three categories of entities:

**Autoantibodies (24 names):**

```python
AUTOANTIBODY_NAMES = [
    "ANA", "anti-dsDNA", "anti-Smith", "RF", "anti-CCP", "anti-Scl-70",
    "anti-centromere", "anti-SSA", "anti-SSB", "anti-Ro", "anti-La",
    "anti-Jo-1", "AChR", "anti-tTG", "TSI", "anti-TPO", "ANCA",
    "anti-cardiolipin", "lupus anticoagulant", "anti-beta2-glycoprotein",
    "anti-RNP", "anti-histone", "anti-Pm-Scl", "anti-RNA Polymerase III",
    "anti-MuSK", "c-ANCA", "p-ANCA", "PR3", "MPO",
]
```

For each detected antibody, the processor extracts:
- **Value**: Numeric value via regex `antibody_name[:\s]*([<>]?\s*\d+\.?\d*)`
- **Positivity**: Boolean via keyword scan (positive, detected, reactive, elevated, abnormal)
- **Titer**: String via regex `antibody_name[^.]*?(1:\d+)`

**Lab tests (45 patterns):**

`LAB_TEST_PATTERNS` is a dictionary mapping test names to regex patterns. Examples:

```python
"CRP": r"c[\s-]?reactive\s+protein|crp\b",
"complement_C3": r"complement\s+c3|c3\s+level|c3\b",
"neurofilament_light": r"neurofilament\s+light|nfl\b|nf[\s-]?l\b",
"NT_proBNP": r"nt[\s-]?pro[\s-]?bnp\b|n[\s-]?terminal[\s-]?pro[\s-]?bnp",
```

Each pattern is compiled with `re.IGNORECASE`. The regex captures the numeric value and optional unit following the test name.

**Dates (4 regex patterns):**

```python
patterns = [
    r"(?:date|visit|collected|drawn)[\s:]*(\d{1,2}[/-]\d{1,2}[/-]\d{2,4})",
    r"(\d{4}-\d{2}-\d{2})",
    r"((?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\w*\s+\d{1,2},?\s+\d{4})",
]
```

Patterns are tried in order. The first match wins. This prioritizes contextual dates ("Visit date: 01/15/2025") over bare ISO dates, which could be document metadata rather than clinical dates.

### 4.4 Semantic Chunking (2500 chars, 200 overlap)

```python
def chunk_text(self, text: str) -> List[str]:
    sentences = re.split(r'(?<=[.!?])\s+', text)
    chunks = []
    current = []
    current_len = 0

    for sentence in sentences:
        if current_len + len(sentence) > self.max_chunk_size and current:
            chunks.append(" ".join(current))
            overlap_text = " ".join(current)
            overlap_start = max(0, len(overlap_text) - self.chunk_overlap)
            overlap = overlap_text[overlap_start:]
            current = [overlap] if overlap else []
            current_len = len(overlap)
        current.append(sentence)
        current_len += len(sentence) + 1
```

The chunker splits on sentence boundaries, then fills chunks up to `max_chunk_size=2500` characters. When a chunk is full, the last 200 characters of the previous chunk carry over as overlap. This ensures that entities or relationships spanning a chunk boundary appear in at least one complete chunk.

### 4.5 Document Classification

**7 doc types** are detected via regex pattern matching:

```python
DOC_TYPE_PATTERNS = {
    "lab_report": ["lab(?:oratory)?\\s+report", "test\\s+results?", "cbc\\b", ...],
    "progress_note": ["progress\\s+note", "office\\s+visit", "chief\\s+complaint", ...],
    "imaging_report": ["radiology", "imaging\\s+report", "mri\\b", ...],
    "pathology_report": ["pathology", "biopsy", "histolog", ...],
    "genetic_report": ["hla\\s+typing", "genetic\\s+test", "pharmacogenomic", ...],
    "referral_letter": ["referral", "dear\\s+(?:dr|doctor)", ...],
    "medication_list": ["medication\\s+(?:list|reconciliation)", "prescription", ...],
}
```

**11 specialties** use a similar pattern-matching approach:

```python
SPECIALTY_PATTERNS = {
    "rheumatology": ["rheumatol", "arthritis", "lupus", "sle\\b", "autoimmun"],
    "neurology": ["neurolog", "ms\\b", "multiple\\s+sclerosis", ...],
    "dermatology": ["dermatol", "skin", "rash", "psoriasis", ...],
    "nephrology": ["nephrol", "kidney", "renal", "proteinuria", ...],
    "gastroenterology": ["gastro", "gi\\b", "bowel", "crohn", ...],
    "ophthalmology": ["ophthal", "eye", "uveitis", "schirmer"],
    "endocrinology": ["endocrin", "thyroid", "diabetes", "graves", ...],
    "cardiology": ["cardiol", "heart", "pots\\b", "tachycardia", ...],
    "allergy_immunology": ["allerg", "immunol", "mast\\s+cell", ...],
    "pulmonology": ["pulmon", "lung", "ild\\b", "pfts?\\b", ...],
    "primary_care": ["primary\\s+care", "pcp\\b", "family\\s+medicine", ...],
}
```

Classification uses a scoring approach: each matching pattern adds 1 point, and the type/specialty with the highest score wins. If no patterns match, the fallback is `"clinical_note"` for doc type and `"general"` for specialty.

### 4.6 Worked Example: Adding a New Document Type

To add "operative_report" as a new document type:

```python
# In document_processor.py, add to DOC_TYPE_PATTERNS:
"operative_report": [
    r"operative\s+report", r"surgical\s+(?:note|report)", r"procedure\s+(?:note|performed)",
    r"anesthesia", r"post[\s-]?operative", r"surgical\s+findings",
],
```

No other changes are needed. The classification system automatically includes the new type in its scoring. Downstream, the `doc_type` field in inserted records will contain `"operative_report"` when matched, and the RAG engine treats all doc types equally during retrieval.

---

## Chapter 5: Extending the Knowledge Base

### 5.1 The Five Knowledge Domains

The knowledge base (`src/knowledge.py`, 855 lines) contains five structured dictionaries:

| Dictionary | Type | Size | Purpose |
|-----------|------|------|---------|
| `HLA_DISEASE_ASSOCIATIONS` | `Dict[str, List[Dict]]` | 22 alleles | HLA allele to disease risk mapping |
| `DISEASE_ACTIVITY_THRESHOLDS` | `Dict[str, Dict]` | 20 scoring systems | Activity score interpretation |
| `AUTOANTIBODY_DISEASE_MAP` | `Dict[str, List[Dict]]` | 24 autoantibodies | Antibody to disease associations |
| `BIOLOGIC_THERAPIES` | `List[Dict]` | 22 therapies | Drug database with PGx |
| `FLARE_BIOMARKER_PATTERNS` | `Dict[str, Dict]` | 13 diseases | Flare prediction biomarker patterns |

### 5.2 Adding a New HLA Association

```python
# In knowledge.py, add to HLA_DISEASE_ASSOCIATIONS:
"HLA-DRB1*16:01": [
    {
        "disease": "systemic_lupus_erythematosus",
        "odds_ratio": 2.0,
        "pmid": "19864127",
        "note": "SLE susceptibility in European populations",
    },
],
```

Required fields: `disease` (must match an `AutoimmuneDisease` enum value or a known disease string), `odds_ratio` (float). Optional fields: `pmid`, `note`.

### 5.3 Adding a New Autoantibody

```python
# In knowledge.py, add to AUTOANTIBODY_DISEASE_MAP:
"anti-SRP": [
    {
        "disease": "necrotizing_autoimmune_myopathy",
        "sensitivity": 0.20,
        "specificity": 0.99,
        "note": "Severe rapidly progressive proximal weakness with necrosis on biopsy",
    },
],
```

Also update `AUTOANTIBODY_NAMES` in `document_processor.py` to enable extraction:

```python
AUTOANTIBODY_NAMES = [
    # ... existing entries ...
    "anti-SRP",
]
```

And update `KNOWLEDGE_VERSION["stats"]["autoantibodies"]` to reflect the new count.

### 5.4 Adding a New Biologic Therapy

```python
# In knowledge.py, add to BIOLOGIC_THERAPIES:
{
    "drug_name": "Anifrolumab",
    "drug_class": "Type I interferon receptor inhibitor",
    "mechanism": "Human anti-IFNAR1 monoclonal antibody -- blocks type I IFN signaling",
    "indicated_diseases": ["systemic_lupus_erythematosus"],
    "pgx_considerations": [
        "IFN gene signature high patients show greater benefit",
        "May affect anti-drug antibody formation",
    ],
    "contraindications": ["Active serious infections", "Active TB"],
    "monitoring_requirements": [
        "Monitor for herpes zoster reactivation",
        "Respiratory tract infection monitoring",
        "SLEDAI-2K and BILAG response assessment",
    ],
},
```

Required fields: `drug_name`, `drug_class`, `indicated_diseases`. The `indicated_diseases` list must use disease strings matching the values in `AutoimmuneDisease` enum.

### 5.5 Adding a New Disease Activity Score

```python
# In knowledge.py, add to DISEASE_ACTIVITY_THRESHOLDS:
"BILAG-2004": {
    "disease": "systemic_lupus_erythematosus",
    "thresholds": {"remission": 0, "low": 1, "moderate": 5, "high": 12},
    "range": [0, 72],
    "components": ["constitutional", "mucocutaneous", "neuropsychiatric",
                   "musculoskeletal", "cardiorespiratory", "renal",
                   "gastrointestinal", "ophthalmic", "haematological"],
    "reference": "PMID:15479896",
},
```

### 5.6 Versioning (v2.0.0 pattern)

The knowledge base is versioned in `KNOWLEDGE_VERSION`:

```python
KNOWLEDGE_VERSION = {
    "version": "2.0.0",
    "last_updated": "2026-03-10",
    "sources": [
        "ACR/EULAR Classification Criteria (2010-2019)",
        "HLA Disease Association Database (PMID:28622507)",
        # ... 7 more sources ...
    ],
    "stats": {
        "hla_alleles": 22,
        "autoantibodies": 24,
        "biologic_therapies": 22,
        "disease_activity_scores": 20,
        "flare_patterns": 13,
        "classification_criteria": 10,
        "overlap_syndromes": 9,
        "lab_test_patterns": 45,
    },
}
```

When modifying the knowledge base:
1. Bump the version (semver: major for breaking changes, minor for additions, patch for corrections).
2. Update `last_updated`.
3. Add any new source references.
4. Update the stats counts.

The version is included in exported reports (Markdown and FHIR) for traceability.

---

## Chapter 6: The Diagnostic Engine

The diagnostic engine (`src/diagnostic_engine.py`, 519 lines) implements clinical decision-support logic separate from the RAG pipeline.

### 6.1 Classification Criteria Evaluation (ACR/EULAR)

The engine supports 10 classification criteria sets, defined in `CLASSIFICATION_CRITERIA`:

| Disease | Criteria Set | Threshold |
|---------|-------------|-----------|
| Rheumatoid Arthritis | 2010 ACR/EULAR RA | >= 6 points |
| Systemic Lupus | 2019 ACR/EULAR SLE | >= 10 points + ANA >= 1:80 |
| Ankylosing Spondylitis | ASAS Axial SpA | >= 1 (imaging or clinical arm) |
| Systemic Sclerosis | 2013 ACR/EULAR SSc | >= 9 points |
| Sjogren's Syndrome | 2016 ACR/EULAR SS | >= 4 points |
| Multiple Sclerosis | 2017 McDonald Criteria | >= 2 (DIS + DIT) |
| Myasthenia Gravis | Clinical Diagnostic Criteria | >= 3 points |
| Celiac Disease | ESPGHAN Criteria | >= 3 points |
| IBD | Montreal Classification | >= 3 points |
| Psoriasis | Clinical Diagnostic Criteria | >= 3 points |

The `evaluate_classification_criteria()` method takes a disease and clinical data dictionary:

```python
def evaluate_classification_criteria(
    self,
    disease: AutoimmuneDisease,
    clinical_data: Dict[str, Any],
) -> Dict[str, Any]:
    criteria_set = CLASSIFICATION_CRITERIA.get(disease)
    total_points = 0
    met_criteria = []
    unmet_criteria = []

    for category, items in criteria_set["criteria"].items():
        if isinstance(items, dict):
            for item, points in items.items():
                if clinical_data.get(item):
                    total_points += points
                    met_criteria.append(f"{category}: {item} (+{points})")
                else:
                    unmet_criteria.append(f"{category}: {item} ({points} pts)")

    meets_criteria = total_points >= criteria_set["threshold"]
```

Usage example:

```python
diag = DiagnosticEngine()
result = diag.evaluate_classification_criteria(
    AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
    {
        "4-10_small": True,           # 3 points (joint involvement)
        "high_positive_RF_or_CCP": True,  # 3 points (serology)
        "abnormal_CRP_or_ESR": True,  # 1 point (acute phase)
        ">=6_weeks": True,            # 1 point (duration)
    }
)
# result["total_points"] = 8, result["meets_criteria"] = True (threshold=6)
```

### 6.2 Differential Diagnosis Algorithm

The `generate_differential()` method scores diseases from two evidence sources:

**Autoantibody scoring:**
```python
for ab in positive_antibodies:
    for assoc in AUTOANTIBODY_DISEASE_MAP.get(ab, []):
        score = assoc.get("specificity", 0.5) * 2.0
        disease_scores[disease] += score
```

Specificity is weighted because highly specific antibodies (anti-CCP for RA at 0.95, anti-dsDNA for SLE at 0.95) are more diagnostically valuable than sensitive but non-specific ones (ANA for SLE at 0.65).

**HLA scoring:**
```python
for allele in hla_alleles:
    for assoc in HLA_DISEASE_ASSOCIATIONS.get(f"HLA-{allele}", []):
        or_score = math.log2(max(assoc["odds_ratio"], 1.0)) * 0.5
        disease_scores[disease] += or_score
```

HLA contributions use log2 of the odds ratio, scaled by 0.5. This logarithmic transform prevents extreme ORs (e.g., HLA-B*27:05 OR=87.4 for AS) from completely dominating the differential. The log2 values:

| Allele | OR | log2(OR) * 0.5 |
|--------|-----|----------------|
| HLA-B*27:05 | 87.4 | 3.22 |
| HLA-C*06:02 | 10.0 | 1.66 |
| HLA-DRB1*04:01 | 4.2 | 1.04 |
| HLA-DRB1*03:01 (SLE) | 2.4 | 0.63 |

### 6.3 Overlap Syndrome Detection

The `detect_overlap_syndromes()` method checks for 9 defined overlap patterns in `OVERLAP_SYNDROMES`:

```python
OVERLAP_SYNDROMES = {
    "mixed_connective_tissue_disease": {
        "required": ["anti-RNP"],
        "features_from": [SLE, SSc, RA],
        "min_features": 2,
    },
    "pots_eds_mcas_triad": {
        "components": ["POTS", "hEDS", "MCAS"],
        "diagnostic_markers": ["tilt_table_positive", "beighton_score_>=5", "tryptase_elevated"],
    },
    # ... 7 more patterns ...
}
```

Detection logic:
1. Check required antibodies (e.g., anti-RNP for MCTD). If any required antibody is absent, skip.
2. Count shared marker matches from the patient's positive antibodies.
3. Assign confidence: "high" if >= 2 markers match, "moderate" if 1 matches.

### 6.4 Adding New Classification Criteria

```python
# In diagnostic_engine.py, add to CLASSIFICATION_CRITERIA:
AutoimmuneDisease.GRAVES_DISEASE: {
    "name": "ATA Guidelines for Graves' Disease Diagnosis",
    "threshold": 3,
    "criteria": {
        "clinical": {
            "hyperthyroidism_symptoms": 1,
            "diffuse_goiter": 1,
            "ophthalmopathy": 2,
        },
        "laboratory": {
            "suppressed_TSH": 2,
            "elevated_free_T4_or_T3": 1,
            "TSI_or_TRAb_positive": 2,
        },
        "imaging": {
            "diffuse_uptake_on_thyroid_scan": 1,
        },
    },
},
```

---

## Chapter 7: The Five Clinical Analysis Engines

The `AutoimmuneAgent` class (`src/agent.py`, 437 lines) orchestrates five clinical analysis engines, executed sequentially in `analyze_patient()`:

### 7.1 Autoantibody Interpretation Engine

```python
def interpret_autoantibodies(self, panel: AutoantibodyPanel) -> List[Dict[str, Any]]:
    findings = []
    for result in panel.results:
        if not result.positive:
            continue
        associations = AUTOANTIBODY_DISEASE_MAP.get(result.antibody, [])
        for assoc in associations:
            findings.append({
                "antibody": result.antibody,
                "disease": assoc["disease"],
                "sensitivity": assoc.get("sensitivity", 0),
                "specificity": assoc.get("specificity", 0),
                "value": result.value,
                "titer": result.titer,
                "pattern": result.pattern,
                "note": assoc.get("note", ""),
            })
    return findings
```

Only positive results are interpreted. Each positive antibody can map to multiple diseases (e.g., ANA maps to SLE, Sjogren's, and SSc). The output includes both the antibody test data (value, titer, pattern) and the association metadata (sensitivity, specificity, clinical notes).

### 7.2 HLA Association Engine

```python
def analyze_hla_associations(self, hla_profile: HLAProfile) -> List[Dict[str, Any]]:
    associations = []
    for allele in hla_profile.all_alleles:
        matches = HLA_DISEASE_ASSOCIATIONS.get(f"HLA-{allele}", [])
        # Also check broader allele groups (e.g., B*27:05 -> B*27)
        if not matches and ":" in allele:
            broad = allele.split(":")[0]
            for key, assocs in HLA_DISEASE_ASSOCIATIONS.items():
                if broad in key:
                    matches.extend(assocs)
    associations.sort(key=lambda x: x["odds_ratio"], reverse=True)
```

The engine performs both exact allele matching (e.g., `HLA-B*27:05`) and broad allele group matching (e.g., `B*27` matching any `B*27:xx` subtype). Results are sorted by odds ratio (highest risk first).

The `HLAProfile` model aggregates alleles from five loci:

```python
@property
def all_alleles(self) -> List[str]:
    return self.hla_a + self.hla_b + self.hla_c + self.hla_drb1 + self.hla_dqb1
```

### 7.3 Disease Activity Scoring Engine

```python
def calculate_disease_activity(self, biomarkers, conditions) -> List[DiseaseActivityScore]:
```

For each diagnosed condition, the engine finds applicable scoring systems from `DISEASE_ACTIVITY_THRESHOLDS` (20 scoring systems across 10 diseases). It uses available biomarkers (CRP, ESR) to estimate activity level:

| Level | Condition |
|-------|-----------|
| REMISSION | marker_value < threshold["remission"] |
| LOW | marker_value < threshold["low"] |
| MODERATE | marker_value < threshold["moderate"] |
| HIGH | marker_value >= threshold["moderate"] |

The engine supports 20 scoring systems: DAS28-CRP, DAS28-ESR, SLEDAI-2K, CDAI, BASDAI, SDAI, PASI, Mayo Score, Harvey-Bradshaw Index, ESSDAI, mRSS, EDSS, QMGS, Marsh Score, Burch-Wartofsky Score, ASDAS, MG-ADL, DAPSA, HbA1c-T1D, and TSH-Hashimoto.

### 7.4 Flare Prediction Engine

The flare prediction algorithm starts with a base risk of 0.3 and adjusts based on biomarker values:

```python
risk_score = 0.3  # Base risk

for marker in pattern.get("early_warning_biomarkers", []):
    value = biomarkers.get(marker, biomarkers.get(marker.lower()))
    if value is not None:
        if marker in ("CRP", "ESR", "IL-6", "calprotectin") and value > 5:
            contributing.append(f"Elevated {marker}: {value}")
            risk_score += 0.15
        elif marker in ("complement_C3", "complement_C4") and value < 80:
            contributing.append(f"Low {marker}: {value}")
            risk_score += 0.15
        elif marker == "albumin" and value < 3.5:
            contributing.append(f"Low albumin: {value}")
            risk_score += 0.1

risk_score = min(max(risk_score, 0.0), 1.0)
```

Risk thresholds (configurable via settings):

| Risk Level | Threshold | Default |
|-----------|-----------|---------|
| IMMINENT | >= `FLARE_RISK_IMMINENT` | >= 0.8 |
| HIGH | >= `FLARE_RISK_HIGH` | >= 0.6 |
| MODERATE | >= `FLARE_RISK_MODERATE` | >= 0.4 |
| LOW | < `FLARE_RISK_MODERATE` | < 0.4 |

Worked example: RA patient with CRP=12 and ESR=45:
- Base risk: 0.3
- Elevated CRP (>5): +0.15 -> 0.45
- Elevated ESR (>5): +0.15 -> 0.60
- Final risk: 0.60 -> FlareRisk.HIGH

### 7.5 Biologic Therapy Recommendation Engine (with PGx filtering)

```python
def recommend_biologics(self, conditions, genotypes=None) -> List[BiologicTherapy]:
    for therapy_data in BIOLOGIC_THERAPIES:
        indicated = therapy_data.get("indicated_diseases", [])
        if not any(c in indicated for c in condition_strs):
            continue
        therapy = BiologicTherapy(
            drug_name=therapy_data["drug_name"],
            drug_class=therapy_data["drug_class"],
            mechanism=therapy_data.get("mechanism", ""),
            pgx_considerations=therapy_data.get("pgx_considerations", []),
            contraindications=therapy_data.get("contraindications", []),
            monitoring_requirements=therapy_data.get("monitoring_requirements", []),
        )
        recommendations.append(therapy)
```

The engine filters the 22 biologic therapies by indication match. PGx considerations are attached to each recommendation (e.g., "CYP3A4 and CYP2C19 metabolism" for tofacitinib, "FCGR3A V158F affects ADCC" for rituximab).

Drug classes represented: TNF inhibitors (5), IL-6R inhibitors (2), Anti-CD20 (2), IL-17A inhibitors (2), IL-23 inhibitors (2), IL-12/23 inhibitor (1), JAK inhibitors (3), BLyS inhibitor (1), T-cell co-stimulation modulator (1), Integrin inhibitors (2), TYK2 inhibitor (1).

---

## Chapter 8: Export System Deep Dive

The export module (`src/export.py`, 389 lines) generates clinical reports in three formats.

### 8.1 Markdown Export

```python
def to_markdown(self, patient_id, analysis_result=None, query_answer=None, evidence_hits=None) -> str:
```

The Markdown report follows a structured template:
1. Header with patient ID, generation timestamp, and knowledge base version.
2. Critical alerts section (if any), formatted as bold list items.
3. Disease activity scores in a Markdown table (Score | Value | Level | Disease).
4. Flare risk predictions with contributing/protective factors and monitoring recommendations.
5. HLA-disease associations in a table (Allele | Disease | Odds Ratio | PMID).
6. Biologic therapy recommendations with mechanism, PGx considerations, contraindications, and monitoring.
7. Evidence sources with relevance badges (green/yellow/red).
8. Footer with clinical validation disclaimer.

### 8.2 FHIR R4 Export

```python
def to_fhir_r4(self, patient_id, analysis_result=None, query_answer=None) -> Dict[str, Any]:
```

The FHIR R4 export produces a Bundle containing:

- **Patient** resource with identifier
- **DiagnosticReport** with:
  - Status: "final"
  - Category: LAB (from v2-0074 code system)
  - Code: 11502-2 (Laboratory report, LOINC)
  - Subject reference to Patient
  - Conclusion summarizing all findings
- **Observation** resources for each disease activity score:
  - Code text: "{score_name} ({disease})"
  - Value: score as quantity
  - Interpretation: activity level
- **Observation** resources for each flare risk prediction:
  - Code text: "Flare Risk Prediction ({disease})"
  - Value: risk score as probability
  - Interpretation: risk level

### 8.3 PDF Export via ReportLab

```python
def to_pdf_bytes(self, patient_id, analysis_result=None, ...) -> bytes:
```

The PDF uses ReportLab with NVIDIA brand colors (`#76B900` for headers). Elements:
- `SimpleDocTemplate` with letter pagesize and 0.75-inch margins
- Custom `ParagraphStyle` instances for title, heading, body, and alert text
- `Table` with alternating row colors for disease activity scores
- Structured paragraphs for biologic therapy recommendations
- Footer with generation metadata

If ReportLab is not installed, the method falls back to returning the Markdown text encoded as UTF-8 bytes.

### 8.4 Adding a New Export Format

To add CSV export:

```python
def to_csv(self, patient_id: str, analysis_result=None) -> str:
    """Export disease activity scores as CSV."""
    import csv
    import io
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["Patient ID", "Score Name", "Value", "Level", "Disease"])
    for score in getattr(analysis_result, "disease_activity_scores", []):
        disease = score.disease.value if hasattr(score.disease, "value") else str(score.disease)
        level = score.level.value if hasattr(score.level, "value") else str(score.level)
        writer.writerow([patient_id, score.score_name, score.score_value, level, disease])
    return output.getvalue()
```

Then add the format to the `/export` endpoint in `api/main.py`:

```python
elif req.format == "csv":
    csv_text = exporter.to_csv(req.patient_id, analysis_result=result)
    return {"format": "csv", "content": csv_text}
```

---

## Chapter 9: Timeline Builder and Diagnostic Odyssey

The timeline builder (`src/timeline_builder.py`, 251 lines) constructs patient diagnostic timelines from ingested clinical documents.

### 9.1 Event Types (12 types)

```python
EVENT_PATTERNS = {
    "symptom_onset":    [r"(?:first|initial|new)\s+(?:complaint|symptom|presentation)", ...],
    "diagnosis":        [r"(?:diagnosed|diagnosis)\s+(?:of|with|:)", ...],
    "misdiagnosis":     [r"(?:previously|initially)\s+(?:diagnosed|labeled|treated)", ...],
    "lab_result":       [r"lab(?:oratory)?\s+result", ...],
    "imaging":          [r"(?:x[\s-]?ray|mri|ct|ultrasound)\s+(?:shows?|reveals?)", ...],
    "biopsy":           [r"biopsy\s+(?:shows?|reveals?|confirms?)", ...],
    "genetic_test":     [r"hla\s+(?:typing|test|result)", ...],
    "treatment_start":  [r"(?:started|initiated|began|prescribed)\s+(?:on\s+)?(?:\w+mab)", ...],
    "treatment_change": [r"(?:switched|changed|transitioned)\s+(?:to|from)", ...],
    "flare":            [r"flare[\s-]?up", r"(?:disease|symptom)\s+exacerbation", ...],
    "referral":         [r"referred?\s+to\s+(?:\w+\s+)?(?:rheumatol|neurolog)", ...],
    "er_visit":         [r"emergency\s+(?:room|department|visit)", ...],
}
```

Each event type has 2-3 regex patterns. Classification scores by counting matches and selecting the highest-scoring type. If no patterns match, the fallback is `"clinical_note"`.

### 9.2 Event Extraction from Documents

```python
def extract_events_from_chunks(self, chunks, patient_id="") -> List[Dict[str, Any]]:
    events = []
    for chunk in chunks:
        text = chunk.get("text_chunk", chunk.get("text", ""))
        event_type = self.classify_event(text)
        event_date = self.extract_date(text) or chunk.get("visit_date", "")
        description = self._summarize_event(text, event_type)
        events.append({...})
```

Each document chunk produces one event. The event date comes from either the text content (via date extraction) or the document metadata (`visit_date` field).

### 9.3 Date Parsing Strategies

Four date patterns are tried in order:

```python
DATE_PATTERNS = [
    (r"(\d{4}-\d{2}-\d{2})", "%Y-%m-%d"),                          # ISO: 2025-01-15
    (r"(\d{1,2}/\d{1,2}/\d{4})", "%m/%d/%Y"),                     # US: 01/15/2025
    (r"(\d{1,2}/\d{1,2}/\d{2})", "%m/%d/%y"),                     # Short: 01/15/25
    (r"((?:Jan|Feb|...)\w*\s+\d{1,2},?\s+\d{4})", "%B %d, %Y"),  # Written: January 15, 2025
]
```

All extracted dates are normalized to ISO format (`%Y-%m-%d`) for consistent sorting.

### 9.4 Temporal Ordering and Pattern Detection

Events are sorted chronologically:

```python
events.sort(key=lambda e: e.get("event_date", "9999"))
```

Events without dates sort to the end (date "9999"). After sorting, the builder assigns `days_from_first_symptom` to each event by finding the earliest `symptom_onset` event and computing deltas.

The `build_timeline()` method aggregates statistics:

```python
return {
    "patient_id": patient_id,
    "total_events": len(events),
    "events": events,
    "specialties_seen": sorted(specialties),
    "event_type_counts": event_types,       # e.g., {"lab_result": 12, "referral": 5, ...}
    "date_range": {"first": ..., "last": ...},
}
```

### 9.5 Misdiagnosis and Delay Detection

The `DiagnosticEngine.analyze_diagnostic_odyssey()` method computes:

- **Diagnostic delay**: Days/months/years from first symptom to diagnosis.
- **Specialist count**: Number of distinct specialties visited.
- **Misdiagnosis list**: Each with date, wrong diagnosis description, and provider.
- **Key diagnostic tests**: Lab results, imaging, biopsies, and genetic tests that contributed to the final diagnosis.

```python
if first_symptom and diagnosis_date:
    d1 = datetime.strptime(first_symptom[:10], "%Y-%m-%d")
    d2 = datetime.strptime(diagnosis_date[:10], "%Y-%m-%d")
    delay_days = (d2 - d1).days
    delay_info = {
        "days": delay_days,
        "months": round(delay_days / 30.44, 1),
        "years": round(delay_days / 365.25, 1),
    }
```

This is clinically significant: the average diagnostic delay for SLE is 4.6 years, and for POTS it is 5-7 years. Quantifying the delay and identifying misdiagnoses provides actionable insights for improving diagnostic pathways.

---

## Chapter 10: Testing Strategies

### 10.1 Test Architecture (8 files, 455 tests)

```
tests/
  test_autoimmune.py          # Core agent tests
  test_api.py                 # FastAPI endpoint tests
  test_collections.py         # Collection manager tests
  test_diagnostic_engine.py   # Diagnostic engine tests
  test_export.py              # Export format tests
  test_rag_engine.py          # RAG engine tests
  test_timeline_builder.py    # Timeline builder tests
  test_production_readiness.py # Production readiness checks
```

### 10.2 Unit Test Patterns

The codebase follows a consistent test pattern:

```python
class TestAutoantibodyInterpretation:
    def test_positive_ana_maps_to_sle(self):
        agent = AutoimmuneAgent()
        panel = AutoantibodyPanel(
            patient_id="test_001",
            collection_date="2026-01-15",
            results=[
                AutoantibodyResult(
                    antibody="ANA",
                    value=320,
                    unit="titer",
                    positive=True,
                    titer="1:320",
                    pattern="homogeneous",
                )
            ],
        )
        findings = agent.interpret_autoantibodies(panel)
        diseases = [f["disease"] for f in findings]
        assert "systemic_lupus_erythematosus" in diseases

    def test_negative_antibody_ignored(self):
        agent = AutoimmuneAgent()
        panel = AutoantibodyPanel(
            patient_id="test_002",
            collection_date="2026-01-15",
            results=[
                AutoantibodyResult(antibody="anti-CCP", value=5.0, positive=False)
            ],
        )
        findings = agent.interpret_autoantibodies(panel)
        assert len(findings) == 0
```

### 10.3 Testing Without Milvus

The RAG engine and collection manager can be tested without a running Milvus instance by mocking:

```python
from unittest.mock import MagicMock, patch

def test_rag_retrieve_without_milvus():
    mock_cm = MagicMock()
    mock_cm.search_all.return_value = {
        "autoimmune_autoantibody_panels": [
            {
                "id": "test_hit_001",
                "score": 0.85,
                "text_chunk": "ANA positive 1:640 homogeneous pattern...",
            }
        ]
    }

    mock_embedder = MagicMock()
    mock_embedder.encode.return_value = MagicMock(tolist=lambda: [0.1] * 384)

    mock_settings = MagicMock()
    mock_settings.TOP_K_PER_COLLECTION = 5
    mock_settings.SCORE_THRESHOLD = 0.40
    mock_settings.MAX_EVIDENCE_ITEMS = 30
    mock_settings.BGE_INSTRUCTION = "Represent this sentence for searching relevant passages: "
    mock_settings.CITATION_HIGH = 0.80
    mock_settings.CITATION_MEDIUM = 0.60
    mock_settings.collection_config = {
        "autoimmune_autoantibody_panels": {"weight": 0.12, "label": "Autoantibody"},
    }

    engine = AutoimmuneRAGEngine(
        collection_manager=mock_cm,
        embedder=mock_embedder,
        llm_client=None,
        settings=mock_settings,
        knowledge=True,
    )

    result = engine.retrieve("What does a positive ANA homogeneous pattern indicate?")
    assert len(result.hits) == 1
    assert result.hits[0].relevance == "high"
```

### 10.4 API Testing with FastAPI TestClient

```python
from fastapi.testclient import TestClient

def test_health_endpoint():
    from api.main import app
    client = TestClient(app)
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "healthy"
    assert "milvus_connected" in data

def test_differential_endpoint():
    client = TestClient(app)
    response = client.post("/differential", json={
        "positive_antibodies": ["anti-CCP", "RF"],
        "hla_alleles": ["DRB1*04:01"],
    })
    assert response.status_code == 200
    data = response.json()
    assert data["differential"][0]["disease"] == "rheumatoid_arthritis"
```

### 10.5 Production Readiness Tests

The `test_production_readiness.py` file validates deployment requirements:

- All required environment variables are documented
- Collection schemas match expected field counts
- Knowledge base version is current
- API endpoints return expected status codes
- Settings validation catches weight sum errors
- CORS configuration is properly restrictive
- API key authentication works when enabled

---

## Chapter 11: Performance Optimization

### 11.1 Parallel Search Tuning (ThreadPoolExecutor across 14 collections)

The `search_all()` method uses `max_workers=6` by default. This parameter is tunable:

```python
def search_all(self, ..., max_workers: int = 6) -> Dict[str, List[Dict]]:
```

Tuning considerations:
- **max_workers=1**: Sequential search. Total time = sum of all collection search times. Use only for debugging.
- **max_workers=6**: Default. Good balance for DGX Spark where Milvus handles concurrent queries well but has lock contention above 8 connections.
- **max_workers=14**: One thread per collection. Maximum parallelism but may overwhelm Milvus with concurrent searches.

With 14 collections at ~3ms per search:
- Sequential: ~42ms total
- Parallel (6 workers): ~9ms total (limited by slowest batch)
- Parallel (14 workers): ~5ms total (limited by slowest collection)

### 11.2 Embedding Caching (256-entry LRU)

The `_embed_cache` in `AutoimmuneRAGEngine` stores up to 256 query embeddings:

```python
self._embed_cache: Dict[str, List[float]] = {}
self._embed_cache_max = 256
```

Cache eviction is FIFO (oldest entry removed when cache is full):

```python
if len(self._embed_cache) >= self._embed_cache_max:
    oldest_key = next(iter(self._embed_cache))
    del self._embed_cache[oldest_key]
```

This relies on Python 3.7+ dictionary insertion order guarantee. The cache key is `text[:512]`, meaning queries that differ only after character 512 will collide. Cache sizing: each 384-float embedding uses ~3KB. 256 entries = ~768KB total -- negligible memory for significant latency savings on repeated queries.

### 11.3 Milvus Index Parameters

```python
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}
SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}
```

Tuning:
- **nlist**: Number of clusters. For collections with N vectors, `nlist = 4 * sqrt(N)` is a good starting point. For 50K vectors: `4 * sqrt(50000) = 894`, so 1024 is appropriate.
- **nprobe**: Number of clusters to search. Higher values improve recall at the cost of latency. For clinical applications where recall matters more than milliseconds, values of 16-32 are recommended.

### 11.4 Score Threshold Tuning (0.40 default)

```python
SCORE_THRESHOLD: float = 0.40
```

The threshold filters out low-relevance hits before they reach the LLM. Tuning:
- **0.30**: More permissive. Include tangentially related evidence. Risk: LLM prompt pollution with irrelevant content.
- **0.40**: Default. Good balance for autoimmune clinical queries.
- **0.50**: More restrictive. Only semantically close hits. Risk: missing relevant evidence on rare conditions.

### 11.5 Deduplication (ID + content hash)

Deduplication in `retrieve()` operates at two levels:

```python
seen_ids = set()
seen_texts = set()

for coll_name, coll_hits in raw_results.items():
    for h in coll_hits:
        # Level 1: ID deduplication
        if h["id"] in seen_ids:
            continue
        seen_ids.add(h["id"])

        # Level 2: Content hash deduplication
        text_hash = hashlib.md5(text[:300].encode()).hexdigest()
        if text_hash in seen_texts:
            continue
        seen_texts.add(text_hash)
```

ID deduplication catches exact duplicates. Content hash deduplication catches near-duplicates where the same text appears in different collections with different IDs. The MD5 hash of the first 300 characters is used because clinical texts often have identical openings (e.g., "Patient presents with...") but diverge later.

---

## Chapter 12: Production Deployment

### 12.1 Docker Multi-Stage Build

A typical Dockerfile for the autoimmune agent:

```dockerfile
# Stage 1: Build dependencies
FROM python:3.11-slim AS builder
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir --prefix=/install -r requirements.txt

# Stage 2: Runtime
FROM python:3.11-slim
WORKDIR /app
COPY --from=builder /install /usr/local
COPY . .

# Pre-download embedding model
RUN python -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('BAAI/bge-small-en-v1.5')"

EXPOSE 8531 8532
HEALTHCHECK --interval=30s --timeout=5s --retries=3 CMD curl -f http://localhost:8532/healthz || exit 1

CMD ["python", "-m", "uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8532"]
```

### 12.2 Docker Compose Configuration

From `docker-compose.dgx-spark.yml`:

```yaml
autoimmune-agent:
  build: ./ai_agent_adds/precision_autoimmune_agent
  ports:
    - "8531:8531"  # Streamlit UI
    - "8532:8532"  # FastAPI API
  environment:
    - AUTO_MILVUS_HOST=milvus
    - AUTO_MILVUS_PORT=19530
    - AUTO_ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
    - AUTO_LOG_LEVEL=INFO
  depends_on:
    - milvus
    - etcd
    - minio
  restart: unless-stopped
  deploy:
    resources:
      limits:
        memory: 4G
```

### 12.3 Environment Variable Management

All configuration is via the `AUTO_` prefix:

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `AUTO_MILVUS_HOST` | str | localhost | Milvus server host |
| `AUTO_MILVUS_PORT` | int | 19530 | Milvus server port |
| `AUTO_ANTHROPIC_API_KEY` | str | "" | Claude API key |
| `AUTO_LLM_MODEL` | str | claude-sonnet-4-6 | LLM model identifier |
| `AUTO_LLM_MAX_TOKENS` | int | 4096 | Max response tokens |
| `AUTO_LLM_TEMPERATURE` | float | 0.2 | LLM temperature |
| `AUTO_TOP_K_PER_COLLECTION` | int | 5 | Results per collection |
| `AUTO_SCORE_THRESHOLD` | float | 0.40 | Minimum cosine similarity |
| `AUTO_MAX_EVIDENCE_ITEMS` | int | 30 | Max evidence items for LLM |
| `AUTO_CONVERSATION_MEMORY_SIZE` | int | 3 | Conversation turns to retain |
| `AUTO_MAX_CHUNK_SIZE` | int | 2500 | Characters per text chunk |
| `AUTO_CHUNK_OVERLAP` | int | 200 | Overlap between chunks |
| `AUTO_CITATION_HIGH` | float | 0.80 | High relevance threshold |
| `AUTO_CITATION_MEDIUM` | float | 0.60 | Medium relevance threshold |
| `AUTO_FLARE_RISK_IMMINENT` | float | 0.8 | Imminent flare threshold |
| `AUTO_FLARE_RISK_HIGH` | float | 0.6 | High flare threshold |
| `AUTO_FLARE_RISK_MODERATE` | float | 0.4 | Moderate flare threshold |
| `AUTO_STREAMING_ENABLED` | bool | True | Enable SSE streaming |
| `AUTO_API_KEY` | str | "" | API auth key (empty=no auth) |
| `AUTO_CORS_ORIGINS` | str | localhost:8080,8531 | Allowed CORS origins |
| `AUTO_MAX_REQUEST_SIZE_MB` | int | 50 | Max PDF upload size |
| `AUTO_REQUEST_TIMEOUT_SECONDS` | int | 60 | Request timeout |
| `AUTO_MILVUS_TIMEOUT_SECONDS` | int | 10 | Milvus query timeout |
| `AUTO_LLM_MAX_RETRIES` | int | 3 | LLM retry count |
| `AUTO_METRICS_ENABLED` | bool | True | Enable Prometheus metrics |
| `AUTO_EMBEDDING_MODEL` | str | BAAI/bge-small-en-v1.5 | Embedding model |
| `AUTO_EMBEDDING_DIM` | int | 384 | Embedding dimensions |
| `AUTO_STREAMLIT_PORT` | int | 8531 | Streamlit UI port |
| `AUTO_API_PORT` | int | 8532 | FastAPI port |

### 12.4 Health Checks

Three health endpoints serve different purposes:

| Endpoint | Method | Purpose | Used By |
|----------|--------|---------|---------|
| `/` | GET | Service identity and version | Browsers, discovery |
| `/health` | GET | Detailed health with stats | Monitoring dashboards |
| `/healthz` | GET | Lightweight liveness probe | Landing page, k8s |

The `/health` endpoint returns:

```json
{
    "status": "healthy",
    "service": "autoimmune-agent",
    "milvus_connected": true,
    "collections": 14,
    "total_vectors": 25000,
    "embedder_loaded": true,
    "llm_available": true,
    "uptime_seconds": 3600
}
```

### 12.5 Monitoring with Prometheus/Grafana

The `/metrics` endpoint exposes Prometheus-compatible metrics:

```
# HELP autoimmune_agent_up Whether the agent is running
# TYPE autoimmune_agent_up gauge
autoimmune_agent_up 1

# HELP autoimmune_collection_vectors Number of vectors per collection
# TYPE autoimmune_collection_vectors gauge
autoimmune_collection_vectors{collection="autoimmune_clinical_documents"} 1500
autoimmune_collection_vectors{collection="autoimmune_patient_labs"} 800
...

# HELP autoimmune_agent_uptime_seconds Agent uptime
# TYPE autoimmune_agent_uptime_seconds gauge
autoimmune_agent_uptime_seconds 3600
```

Configure Prometheus scraping in `prometheus.yml`:

```yaml
- job_name: 'autoimmune-agent'
  static_configs:
    - targets: ['localhost:8532']
  metrics_path: '/metrics'
  scrape_interval: 30s
```

---

## Chapter 13: Integration with HCLS AI Factory

### 13.1 The 3-Stage Pipeline

The Precision Autoimmune Agent integrates with the HCLS AI Factory's three-stage precision medicine pipeline:

1. **Genomics Pipeline** (`genomics-pipeline/`): FASTQ -> BAM -> VCF using Parabricks/DeepVariant/BWA-MEM2. Produces variant calls that feed into HLA typing and pharmacogenomic analysis.
2. **RAG/Chat Pipeline** (`rag-chat-pipeline/`): Milvus + Claude AI for variant interpretation. Shares the Milvus infrastructure with the autoimmune agent.
3. **Drug Discovery Pipeline** (`drug-discovery-pipeline/`): BioNeMo MolMIM/DiffDock/RDKit for drug candidate generation. The autoimmune agent's biologic therapy recommendations can guide target selection.

### 13.2 The Genomic Evidence Bridge (shared collection, weight 0.02)

The `genomic_evidence` collection is shared and read-only, populated by the genomics pipeline:
```python
COLL_GENOMIC_EVIDENCE: str = "genomic_evidence"  # shared read-only
WEIGHT_GENOMIC_EVIDENCE: float = 0.02
```

The autoimmune agent reads from this collection but never writes to it. In `create_all_collections()`:

```python
# Skip genomic_evidence — shared, not ours to create
if name == "genomic_evidence":
    if utility.has_collection(name, using=self._alias):
        collections[name] = Collection(name, using=self._alias)
        collections[name].load()
    continue
```

The low weight (0.02) reflects that genomic evidence provides supporting context but is not the primary data source for autoimmune analysis. A patient's VCF-derived variants might identify HLA alleles or pharmacogenomic variants relevant to therapy selection.

### 13.3 Cross-Agent Communication

The `AutoimmuneAgent` class includes three cross-agent integration points:

```python
def request_biomarker_context(self, patient_id, biomarker_names) -> Dict:
    """Request inflammation context from the Biomarker Agent."""
    # Returns trends for CRP, ESR, IL-6, etc.

def request_imaging_context(self, patient_id, body_regions) -> Dict:
    """Request imaging findings from the Imaging Agent."""
    # Returns joint damage scores, organ involvement

def publish_diagnosis_event(self, patient_id, disease, confidence, evidence) -> Dict:
    """Publish diagnosis for other agents to consume."""
    # Event-driven notification to Biomarker, Imaging, Oncology agents
```

Currently these return stub responses. In production, they would communicate via:
- HTTP calls to other agent APIs (Biomarker at port 8530, Imaging TBD)
- Event bus (Redis Streams or Kafka) for asynchronous notification

### 13.4 Event Bus Architecture

The `publish_diagnosis_event()` method emits structured events:

```python
event = {
    "event_type": "autoimmune_diagnosis",
    "source_agent": "precision_autoimmune",
    "patient_id": patient_id,
    "disease": disease,
    "confidence": confidence,
    "supporting_evidence": supporting_evidence,
    "timestamp": None,  # Set by event bus
}
```

When implemented, the event bus enables:
- **Biomarker Agent**: Adjusts monitoring panels based on new autoimmune diagnosis
- **Imaging Agent**: Prioritizes imaging modalities for affected organs
- **Oncology Agent**: Flags autoimmune/immunotherapy contraindications

### 13.5 Shared Milvus Infrastructure

All five intelligence agents share a single Milvus instance (default port 19530) with etcd for metadata and MinIO for object storage:

```yaml
# From docker-compose.dgx-spark.yml
milvus:
  image: milvusdb/milvus:latest
  ports:
    - "19530:19530"
  depends_on:
    - etcd
    - minio
```

Each agent uses its own connection alias to avoid interference:

```python
self._alias = "autoimmune_agent"
connections.connect(alias=self._alias, host=self.host, port=self.port)
```

Collection names are prefixed by domain (`autoimmune_`, `cart_`, `biomarker_`, etc.) to prevent name collisions. The only shared collection is `genomic_evidence`.

---

## Chapter 14: Future Architecture

### 14.1 Multi-Agent Orchestration

The current cross-agent integration uses point-to-point stubs. Future architecture would introduce an orchestrator that coordinates multi-agent analysis:

```
Patient Genome → Genomics Pipeline → Variant Calls
                                          ↓
                                    Orchestrator
                                   /     |       \
                         Biomarker   Autoimmune   Oncology
                          Agent       Agent        Agent
                                   \     |       /
                                    Unified Report
```

The orchestrator would manage query routing (which agent handles which question), result aggregation, and conflict resolution (e.g., when autoimmune and oncology agents recommend conflicting therapies).

### 14.2 Graph Databases for Knowledge

The current knowledge base is a set of Python dictionaries. A future enhancement would migrate to a graph database (Neo4j or Amazon Neptune) to enable:
- **Relationship queries**: "What diseases share HLA-DRB1*03:01 and anti-SSA?"
- **Path finding**: "What is the shortest diagnostic path from ANA+ to confirmed SLE?"
- **Temporal reasoning**: "How do biomarker trajectories predict disease transition from UCTD to SLE?"

### 14.3 Fine-Tuned Domain Embeddings

BGE-small-en-v1.5 is a general-purpose embedding model. Fine-tuning on autoimmune clinical text would improve retrieval:
- Training data: Autoantibody reports, HLA typing results, rheumatology progress notes
- Expected improvement: 5-15% recall improvement on domain-specific queries
- Trade-off: Fine-tuned models must be maintained alongside general models

### 14.4 Real-Time Data Streaming

Current ingestion is batch-oriented (upload PDFs, process, embed, insert). Real-time streaming would enable:
- HL7 FHIR resource streaming from EHR systems
- Continuous lab result ingestion and flare risk re-calculation
- Real-time alert generation when biomarker patterns trigger threshold crossings

### 14.5 VAST AI OS Integration

The HCLS AI Factory is being prepared for deployment on VAST AI OS, which provides:
- Distributed GPU infrastructure for running multiple agents concurrently
- Shared storage for large model weights and embedding indices
- Auto-scaling based on query load
- Deployment templates in `aios/` directory

---

## Appendix A: Complete API Reference

All 14 endpoints exposed by the FastAPI server (`api/main.py`):

| # | Method | Path | Description | Auth |
|---|--------|------|-------------|------|
| 1 | GET | `/` | Service identity, version, ports | No |
| 2 | GET | `/health` | Detailed health check (Milvus, embedder, LLM, uptime) | No |
| 3 | GET | `/healthz` | Lightweight liveness probe | No |
| 4 | GET | `/metrics` | Prometheus-compatible metrics | No |
| 5 | POST | `/query` | Full RAG query (retrieve + synthesize) | Yes |
| 6 | POST | `/query/stream` | Streaming RAG via SSE | Yes |
| 7 | POST | `/search` | Evidence-only search (no LLM) | Yes |
| 8 | POST | `/analyze` | Full patient analysis pipeline | Yes |
| 9 | POST | `/differential` | Differential diagnosis from antibodies/HLA | Yes |
| 10 | POST | `/ingest/upload` | Upload and ingest a clinical PDF | Yes |
| 11 | POST | `/ingest/demo-data` | Ingest all demo patient data | Yes |
| 12 | GET | `/collections` | List collections with vector counts | Yes |
| 13 | POST | `/collections/create` | Create all collections (optional drop) | Yes |
| 14 | POST | `/export` | Export report (markdown, fhir, pdf) | Yes |

**Request/Response examples:**

```bash
# Full RAG query
curl -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your_key" \
  -d '{
    "question": "What are the HLA associations for ankylosing spondylitis?",
    "patient_id": "sarah_mitchell",
    "top_k": 5
  }'

# Response:
{
  "answer": "Ankylosing spondylitis has the strongest known HLA association...",
  "evidence_count": 12,
  "collections_searched": 14,
  "search_time_ms": 45.2
}

# Differential diagnosis
curl -X POST http://localhost:8532/differential \
  -H "Content-Type: application/json" \
  -d '{
    "positive_antibodies": ["ANA", "anti-dsDNA", "anti-Smith"],
    "hla_alleles": ["DRB1*03:01"]
  }'

# Response:
{
  "differential": [
    {"disease": "systemic_lupus_erythematosus", "score": 5.23, "rank": 1, ...},
    {"disease": "sjogrens_syndrome", "score": 1.85, "rank": 2, ...}
  ]
}

# Export to FHIR
curl -X POST http://localhost:8532/export \
  -H "Content-Type: application/json" \
  -d '{"patient_id": "sarah_mitchell", "format": "fhir"}'
```

---

## Appendix B: Configuration Reference

All `AUTO_*` environment variables with types, defaults, and descriptions:

```python
class AutoimmuneSettings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="AUTO_")

    # Paths
    PROJECT_ROOT: Path                          # Auto-detected from settings.py location

    # Milvus
    MILVUS_HOST: str = "localhost"              # Milvus server hostname
    MILVUS_PORT: int = 19530                    # Milvus server port

    # 14 collection name constants
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
    COLL_GENOMIC_EVIDENCE: str = "genomic_evidence"

    # Embedding
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIM: int = 384
    EMBEDDING_BATCH_SIZE: int = 32
    BGE_INSTRUCTION: str = "Represent this sentence for searching relevant passages: "

    # LLM
    ANTHROPIC_API_KEY: str = ""
    LLM_MODEL: str = "claude-sonnet-4-6"
    LLM_MAX_TOKENS: int = 4096
    LLM_TEMPERATURE: float = 0.2

    # RAG parameters
    TOP_K_PER_COLLECTION: int = 5
    SCORE_THRESHOLD: float = 0.40
    MAX_EVIDENCE_ITEMS: int = 30
    CONVERSATION_MEMORY_SIZE: int = 3

    # 14 collection weights (sum ~ 1.0)
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

    # Ports
    STREAMLIT_PORT: int = 8531
    API_PORT: int = 8532

    # Authentication
    API_KEY: str = ""
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8531"
    MAX_REQUEST_SIZE_MB: int = 50

    # Document processing
    MAX_CHUNK_SIZE: int = 2500
    CHUNK_OVERLAP: int = 200
    PDF_DPI: int = 200

    # Relevance thresholds
    CITATION_HIGH: float = 0.80
    CITATION_MEDIUM: float = 0.60

    # Flare risk thresholds
    FLARE_RISK_IMMINENT: float = 0.8
    FLARE_RISK_HIGH: float = 0.6
    FLARE_RISK_MODERATE: float = 0.4

    # Evidence display
    MAX_EVIDENCE_TEXT_LENGTH: int = 1500
    MAX_KNOWLEDGE_CONTEXT_ITEMS: int = 25

    # Streaming
    STREAMING_ENABLED: bool = True

    # Timeouts
    REQUEST_TIMEOUT_SECONDS: int = 60
    MILVUS_TIMEOUT_SECONDS: int = 10
    LLM_MAX_RETRIES: int = 3

    # Logging
    LOG_LEVEL: str = "INFO"
    LOG_DIR: str = ""

    # Metrics
    METRICS_ENABLED: bool = True
```

---

## Appendix C: Metric Reference

Prometheus metrics exposed at `/metrics`:

| Metric | Type | Description |
|--------|------|-------------|
| `autoimmune_agent_up` | gauge | Whether the agent is running (always 1 when reachable) |
| `autoimmune_collection_vectors{collection="..."}` | gauge | Number of vectors in each collection (14 labels) |
| `autoimmune_agent_uptime_seconds` | gauge | Agent uptime in seconds since startup |

Additional metrics can be added by extending the `/metrics` endpoint in `api/main.py`. For richer instrumentation, integrate the `prometheus_client` library with `Histogram` and `Counter` objects for query latency and throughput tracking.

---

*This guide covers the internals of the Precision Autoimmune Intelligence Agent as of version 2.0.0 (March 2026). For the foundations guide, see DEMO_GUIDE.md. For API documentation, see API_REFERENCE.md. For architectural overview, see ARCHITECTURE_GUIDE.md.*

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Autoimmune Agent is a clinical decision support research tool for autoimmune disease assessment. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
