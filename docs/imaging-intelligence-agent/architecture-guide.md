# Imaging Intelligence Agent -- Architecture Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0

---

## 1. System Architecture Overview

The Imaging Intelligence Agent is organized into six interconnected layers, each with clear responsibilities and interfaces. The system is designed for deployment on a single NVIDIA DGX Spark with 128 GB unified memory, but runs equally well in CPU-only mode with mock NIM fallbacks.

### Design Principles

1. **Graceful degradation** -- Every NIM-dependent feature falls back to mock mode automatically
2. **Shared infrastructure** -- Reuses Milvus and embedding model from the HCLS AI Factory platform
3. **Cross-agent interoperability** -- Read-only access to `genomic_evidence` collection from Stage 2
4. **Consistent patterns** -- Follows the same Pydantic BaseSettings, collection manager, and RAG engine patterns as the CAR-T Intelligence Agent

---

## 2. Component Diagram

```
+=========================================================================+
|                        PRESENTATION LAYER                                |
|                                                                          |
|  +---------------------------+    +-------------------------------+      |
|  | Streamlit Chat UI (8525)  |    | FastAPI REST Server (8524)    |      |
|  | app/imaging_ui.py         |    | api/main.py                  |      |
|  |                           |    |   +-- routes/meta_agent.py    |      |
|  | - Chat interface          |    |   +-- routes/nim.py          |      |
|  | - Evidence panel          |    |   +-- routes/workflows.py    |      |
|  | - Workflow runner         |    |   +-- routes/reports.py      |      |
|  | - NIM status              |    |                               |      |
|  | - Report export           |    | Prometheus metrics            |      |
|  +---------------------------+    +-------------------------------+      |
+=========================================================================+
                    |                            |
                    v                            v
+=========================================================================+
|                        INTELLIGENCE LAYER                                |
|                                                                          |
|  +----------------------------+   +-----------------------------+        |
|  | Imaging Agent              |   | RAG Engine                  |        |
|  | src/agent.py               |   | src/rag_engine.py           |        |
|  |                            |   |                             |        |
|  | - Query classification     |   | - Multi-collection search   |        |
|  | - Workflow dispatch        |   | - Weighted scoring          |        |
|  | - NIM orchestration        |   | - Query expansion           |        |
|  | - Response assembly        |   | - Comparative analysis      |        |
|  +----------------------------+   | - LLM synthesis             |        |
|                                   +-----------------------------+        |
|  +----------------------------+   +-----------------------------+        |
|  | Knowledge Graph            |   | Query Expansion             |        |
|  | src/knowledge.py           |   | src/query_expansion.py      |        |
|  |                            |   |                             |        |
|  | - 25 pathologies           |   | - Domain-specific maps      |        |
|  | - 9 modalities             |   | - Keyword -> term expansion |        |
|  | - 21 anatomy entries       |   | - Entity resolution         |        |
|  +----------------------------+   +-----------------------------+        |
+=========================================================================+
                    |                            |
                    v                            v
+=========================================================================+
|                        INFERENCE LAYER                                   |
|                                                                          |
|  +------------------------------------------------------------------+   |
|  | NIM Service Manager (src/nim/service_manager.py)                  |   |
|  |                                                                    |   |
|  |  +-----------+  +-----------+  +-----------+  +----------------+  |   |
|  |  | VISTA-3D  |  | MAISI     |  | VILA-M3   |  | Llama-3 8B     |  |   |
|  |  | Client    |  | Client    |  | Client    |  | Client         |  |   |
|  |  | 8530      |  | 8531      |  | 8532      |  | 8520           |  |   |
|  |  +-----------+  +-----------+  +-----------+  +----------------+  |   |
|  |                                                                    |   |
|  |  All inherit BaseNIMClient: health check + retry + mock fallback  |   |
|  +------------------------------------------------------------------+   |
|                                                                          |
|  +------------------------------------------------------------------+   |
|  | Workflow Engine (src/workflows/)                                   |   |
|  |                                                                    |   |
|  |  +-- CTHeadHemorrhageWorkflow    (< 90 sec, 3D U-Net)            |   |
|  |  +-- CTChestLungNoduleWorkflow   (< 5 min, RetinaNet+SegResNet)  |   |
|  |  +-- CXRRapidFindingsWorkflow    (< 30 sec, DenseNet-121)        |   |
|  |  +-- MRIBrainMSLesionWorkflow    (< 5 min, 3D U-Net+SyN)        |   |
|  |                                                                    |   |
|  |  All inherit BaseImagingWorkflow: preprocess->infer->postprocess  |   |
|  +------------------------------------------------------------------+   |
+=========================================================================+
                    |                            |
                    v                            v
+=========================================================================+
|                        DATA LAYER                                        |
|                                                                          |
|  +----------------------------+   +-----------------------------+        |
|  | Milvus 2.4 (19530)        |   | Collection Manager          |        |
|  |                            |   | src/collections.py          |        |
|  | 10 imaging collections     |   |                             |        |
|  | + 1 read-only genomic      |   | - Schema definitions        |        |
|  | IVF_FLAT / COSINE / 384d   |   | - CRUD operations           |        |
|  +----------------------------+   | - Parallel search            |        |
|                                   +-----------------------------+        |
|  +----------------------------+   +-----------------------------+        |
|  | Pydantic Models            |   | Ingest Pipelines            |        |
|  | src/models.py              |   | src/ingest/                 |        |
|  |                            |   |                             |        |
|  | - 10 collection models     |   | - PubMed parser             |        |
|  | - 4 NIM result models      |   | - ClinicalTrials parser     |        |
|  | - Search result models     |   | - 6 seed data parsers       |        |
|  | - Agent I/O models         |   | - APScheduler integration   |        |
|  +----------------------------+   +-----------------------------+        |
+=========================================================================+
```

---

## 3. Data Flow

### 3.1 RAG Query Flow

```
User Query: "What is ACR Lung-RADS classification?"
       |
       v
[1. Query Classification]
       |-- Detect comparative ("X vs Y")? --> No
       |-- Detect modality filter? --> CT
       |-- Detect body region? --> Chest
       |
       v
[2. Query Expansion]
       |-- "Lung-RADS" --> ["lung_rads", "lung_cancer_screening",
       |                     "nodule_management", "ACR", ...]
       |
       v
[3. Embedding]
       |-- BGE-small-en-v1.5: "Represent this sentence: ..."
       |-- Output: 384-dim float32 vector
       |
       v
[4. Parallel Multi-Collection Search]
       |-- imaging_literature    (weight 0.18, top-5) --> 5 hits
       |-- imaging_guidelines    (weight 0.10, top-5) --> 5 hits
       |-- imaging_findings      (weight 0.15, top-5) --> 3 hits
       |-- imaging_trials        (weight 0.12, top-5) --> 4 hits
       |-- ... (all 11 collections)
       |
       v
[5. Weighted Score Merge]
       |-- Combine hits across collections
       |-- Apply collection weights
       |-- Filter by SCORE_THRESHOLD (0.4)
       |-- Sort by weighted score descending
       |
       v
[6. Knowledge Graph Augmentation]
       |-- Match "lung_nodule" pathology entry
       |-- Inject: Lung-RADS categories, severity criteria,
       |           CT characteristics, AI models
       |
       v
[7. LLM Synthesis]
       |-- Build prompt: question + evidence + knowledge context
       |-- Inject conversation history (up to 3 prior turns)
       |-- Call Claude API (or Llama-3 NIM fallback)
       |
       v
[8. Response Assembly]
       |-- Grounded answer with evidence citations
       |-- Source references with scores
       |-- Follow-up question suggestions
       |-- NIM service availability status
```

### 3.2 Workflow Execution Flow

```
API Request: POST /workflow/ct_head_hemorrhage/run
       |
       v
[1. Workflow Registry Lookup]
       |-- WORKFLOW_REGISTRY["ct_head_hemorrhage"]
       |-- Instantiate CTHeadHemorrhageWorkflow(mock_mode=True)
       |
       v
[2. Preprocess]
       |-- Mock: skip (return synthetic volume metadata)
       |-- Live: LoadImaged -> EnsureChannelFirst -> Orientationd(RAS)
       |         -> Spacingd(1mm) -> ScaleIntensityRanged(0-80 HU)
       |
       v
[3. Infer]
       |-- Mock: return synthetic segmentation result
       |-- Live: 3D U-Net binary segmentation via MONAI
       |
       v
[4. Postprocess]
       |-- Volume estimation: voxel count x voxel volume
       |-- Midline shift: center of mass vs falx cerebri
       |-- Max thickness measurement
       |-- BTF urgency classification (P1/P2/P4)
       |
       v
[5. WorkflowResult]
       |-- findings: [{category, description, severity, recommendation}]
       |-- measurements: {volume_ml, shift_mm, thickness_mm}
       |-- classification: "P1" / "P2" / "P4"
       |-- severity: critical / urgent / routine
       |-- inference_time_ms, is_mock
```

---

## 4. Milvus Collection Design

### 4.1 Index Configuration

All collections use the same index configuration:

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric type | COSINE |
| nlist | 1024 |
| nprobe | 16 |
| Vector dimension | 384 |
| Embedding model | BAAI/bge-small-en-v1.5 |

### 4.2 Schema Pattern

Every collection follows the same field pattern:

```python
FieldSchema("id",        DataType.VARCHAR, max_length=100, is_primary=True)
FieldSchema("embedding", DataType.FLOAT_VECTOR, dim=384)
FieldSchema("text",      DataType.VARCHAR, max_length=3000)
# ... domain-specific metadata fields (VARCHAR, INT64, FLOAT, etc.)
```

### 4.3 Search Strategy

1. **Parallel search:** All collections are searched simultaneously using `ThreadPoolExecutor`
2. **Per-collection top-K:** Default 5 results per collection (configurable 1-50)
3. **Weighted scoring:** Each collection has a configurable weight (0.04 to 0.18)
4. **Score threshold:** Results below 0.4 cosine similarity are filtered out
5. **Asymmetric embedding:** Queries use BGE instruction prefix `"Represent this sentence for searching relevant passages: "`

---

## 5. NIM Client Layer Design

### 5.1 BaseNIMClient (ABC)

```
BaseNIMClient
    |
    +-- health_check()       Ping /v1/health/ready
    +-- is_available()       Cached check (30s interval)
    +-- _request()           HTTP POST with tenacity retry (3 attempts)
    +-- _mock_response()     Abstract: return synthetic result
    +-- _invoke_or_mock()    Try real NIM, fall back to mock
    +-- get_status()         Return "available" / "mock" / "unavailable"
```

### 5.2 Client Hierarchy

```
BaseNIMClient (ABC)
    |
    +-- VISTA3DClient
    |       segment(input_path, target_classes) -> SegmentationResult
    |
    +-- MAISIClient
    |       generate(body_region, resolution) -> SyntheticCTResult
    |
    +-- VILAM3Client
    |       analyze(question, input_path) -> VLMResponse
    |
    +-- LlamaLLMClient
            complete(messages) -> str
            (OpenAI-compatible /v1/chat/completions)
```

### 5.3 NIMServiceManager

Coordinates all four NIM clients:

```python
NIMServiceManager(settings)
    .vista3d    -> VISTA3DClient
    .maisi      -> MAISIClient
    .vilam3     -> VILAM3Client
    .llm        -> LlamaLLMClient
    .check_all_services() -> Dict[str, str]  # name -> status
```

### 5.4 Mock Fallback Logic

```
_invoke_or_mock(endpoint, payload):
    if is_available():
        try:
            return _request(endpoint, payload)    # Real NIM
        except:
            if mock_enabled:
                return _mock_response()           # Fallback mock
            raise
    elif mock_enabled:
        return _mock_response()                   # Direct mock
    else:
        raise ConnectionError
```

---

## 6. Workflow Pipeline Design

### 6.1 BaseImagingWorkflow (ABC)

All four reference workflows inherit from the same abstract base class:

```python
class BaseImagingWorkflow(ABC):
    WORKFLOW_NAME: str
    TARGET_LATENCY_SEC: float
    MODALITY: str
    BODY_REGION: str
    MODELS_USED: List[str]

    preprocess(input_path)    -> Any          # Abstract
    infer(preprocessed)       -> Dict         # Abstract
    postprocess(result)       -> WorkflowResult   # Abstract
    _mock_inference()         -> Dict         # Abstract
    run(input_path)           -> WorkflowResult   # Orchestrator
    get_workflow_info()       -> Dict         # Metadata
```

### 6.2 Workflow Registry

```python
WORKFLOW_REGISTRY = {
    "ct_head_hemorrhage":   CTHeadHemorrhageWorkflow,
    "ct_chest_lung_nodule": CTChestLungNoduleWorkflow,
    "cxr_rapid_findings":   CXRRapidFindingsWorkflow,
    "mri_brain_ms_lesion":  MRIBrainMSLesionWorkflow,
}
```

Dynamic dispatch via the `/workflow/{name}/run` API endpoint.

### 6.3 Error Handling

```
run(input_path):
    start = time.time()
    try:
        if mock_mode:
            raw = _mock_inference()
        else:
            preprocessed = preprocess(input_path)
            raw = infer(preprocessed)
        result = postprocess(raw)
        result.inference_time_ms = elapsed
        result.is_mock = mock_mode
        return result
    except:
        return WorkflowResult(status=FAILED, inference_time_ms=elapsed)
```

---

## 7. Ingest Pipeline Design

### 7.1 Pipeline Pattern

```
[Source] --> fetch() --> parse() --> embed() --> store()
              |            |           |           |
         HTTP/API    Extract fields  BGE-small   Milvus
         PubMed      Normalize       384-dim     upsert
         CT.gov      Validate
         Seed JSON   Pydantic model
```

### 7.2 Ingest Parsers

| Parser | Source | Collection |
|---|---|---|
| `literature_parser.py` | PubMed (NCBI E-utilities) | `imaging_literature` |
| `clinical_trials_parser.py` | ClinicalTrials.gov API v2 | `imaging_trials` |
| `finding_parser.py` | Curated seed data | `imaging_findings` |
| `protocol_parser.py` | Curated seed data | `imaging_protocols` |
| `device_parser.py` | Curated seed data | `imaging_devices` |
| `anatomy_parser.py` | Curated seed data | `imaging_anatomy` |
| `benchmark_parser.py` | Curated seed data | `imaging_benchmarks` |
| `guideline_parser.py` | Curated seed data | `imaging_guidelines` |
| `report_template_parser.py` | Curated seed data | `imaging_report_templates` |

### 7.3 PubMed Client (`src/utils/pubmed_client.py`)

- NCBI E-utilities: esearch + efetch
- Optional API key for increased rate limits
- Configurable max results (default 5000)

### 7.4 Scheduling

APScheduler (`src/scheduler.py`) supports periodic re-ingestion:
- Default interval: 168 hours (weekly)
- Configurable via `IMAGING_INGEST_SCHEDULE_HOURS`
- Disabled by default (`IMAGING_INGEST_ENABLED=false`)

---

## 8. API Layer

### 8.1 FastAPI Application

- **Lifespan management:** Initializes Milvus connection, embedding model, NIM service manager, and RAG engine on startup
- **CORS:** Enabled for all origins (development mode)
- **Prometheus metrics:** Query count, latency histogram, search hit histogram
- **Health check:** Reports collection stats, NIM service status, and overall system health

### 8.2 Route Organization

| Router | Prefix | Tags | Endpoints |
|---|---|---|---|
| `meta_agent` | `/api` | Meta-Agent | `/api/ask` |
| `nim` | `/nim` | NIM Services | `/nim/status`, `/nim/vista3d/segment`, `/nim/maisi/generate`, `/nim/vilam3/analyze` |
| `workflows` | (root) | Workflows | `/workflows`, `/workflow/{name}/info`, `/workflow/{name}/run` |
| `reports` | (root) | Reports | `/reports/generate` |
| `events` | `/events` | DICOM Events | `/events/dicom-webhook`, `/events/history`, `/events/status` |

Core endpoints registered directly on the app: `/health`, `/collections`, `/query`, `/search`, `/find-related`, `/knowledge/stats`, `/metrics`

---

## 9. UI Layer

### 9.1 Streamlit Application (`app/imaging_ui.py`)

The Streamlit UI provides:

1. **Chat interface** with multi-turn conversation memory
2. **Evidence panel** with expandable results grouped by collection
3. **Comparative analysis** auto-detection and dual-panel display
4. **Workflow runner** sidebar for executing reference workflows
5. **NIM service status** indicators showing available/mock/unavailable
6. **Report export** button for PDF generation
7. **Collection statistics** in the sidebar
8. **NVIDIA-themed** dark/green styling

---

## 10. Cross-Modal Integration

### 10.1 Genomic Pipeline Trigger (Implemented)

The `CrossModalTrigger` class (`src/cross_modal.py`) automatically enriches high-risk imaging findings with genomic context from the shared `genomic_evidence` collection (3.5M vectors).

**Trigger conditions:**
- Lung-RADS 4A+ findings --> queries EGFR, ALK, ROS1, KRAS variants
- CXR urgent consolidation --> queries infection-related genomic variants

**Data flow:**
```
WorkflowResult (severity=urgent, classification=Lung-RADS 4A)
    |
    v
CrossModalTrigger.evaluate(workflow_result)
    |
    v
Query genomic_evidence collection (3 queries: EGFR, ALK, KRAS)
    |
    v
CrossModalResult (12 genomic hits, top score: 0.78)
    |
    v
AgentResponse.cross_modal (enriched response)
```

**Configuration:**
```python
CROSS_MODAL_ENABLED: bool = True  # Active
```

### 10.2 FHIR R4 Export Architecture

The `export_fhir()` function generates FHIR R4 DiagnosticReport Bundles:

```
FHIR Bundle (type: collection)
├── Patient resource (stub with identifier)
├── ImagingStudy resource (modality auto-detected from query)
├── Observation resources (one per workflow finding)
│   ├── SNOMED CT coding (finding category)
│   ├── Interpretation (severity → HH/H/A/N)
│   └── Components (measurements with UCUM units)
└── DiagnosticReport resource
    ├── LOINC category (LP29684-5 Radiology)
    ├── LOINC code (18748-4 Diagnostic imaging study)
    ├── conclusionCode (SNOMED for all findings)
    └── extension (cross-modal enrichment summary)
```

### 10.3 DICOM Ingestion Architecture

```
Orthanc DICOM Server (port 8042 HTTP, 4242 C-STORE)
    |
    v
POST /events/dicom-webhook (study.complete event)
    |
    v
determine_workflow(modality, body_region) → workflow name
    |
    v
WorkflowRegistry.run(workflow_name, study_data)
    |
    v
DicomIngestionResult (findings, classification, severity)
    |
    v
Event history (in-memory, max 200 entries)
```

### 10.4 Drug Discovery Pipeline Feed (Phase 2)

```
Quantitative imaging endpoint
    |-- Tumor volume change
    |-- RECIST measurements
    |-- Treatment response
    |
    v
Drug Discovery Pipeline
    |-- Treatment-response tracking
    |-- Molecular target validation
```

---

*For NIM-specific setup instructions, see `NIM_INTEGRATION_GUIDE.md`. For the complete implementation specification, see `PROJECT_BIBLE.md`.*
