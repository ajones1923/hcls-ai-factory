# Clinical Imaging Engine -- Architecture Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 2.1.0

---

## 1. System Architecture Overview

The Clinical Imaging Engine (Engine 4) is organized into six interconnected layers, each with clear responsibilities and interfaces. The system integrates 20 NVIDIA technologies (Community Edition, all free), 9 NIM clients, 9 clinical workflows, 13 Milvus collections (38,028 vectors including 1,938 real PubMed papers), and 1,324 tests. It is designed for deployment on a single NVIDIA DGX Spark ($4,699) with 128 GB unified memory in a 3-tier model (Community/Enterprise/Research), but runs equally well in CPU-only mode with mock NIM fallbacks.

### Design Principles

1. **Graceful degradation** -- Every NIM-dependent feature falls back to mock mode automatically
2. **Shared infrastructure** -- Reuses Milvus and embedding model from the HCLS AI Factory platform
3. **Cross-agent interoperability** -- Read-only access to `genomic_evidence` collection from Stage 2
4. **Consistent patterns** -- Follows the same Pydantic BaseSettings, collection manager, and RAG engine patterns as the CAR-T Intelligence Engine

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
|  |  | VISTA-3D  |  | MAISI     |  | VILA-M3   |  | LLM            |  |   |
|  |  | Client    |  | Client    |  | Client    |  | Client         |  |   |
|  |  | 8530      |  | 8531      |  | 8532      |  | 8520           |  |   |
|  |  +-----------+  +-----------+  +-----------+  +----------------+  |   |
|  |  +-----------+  +-----------+  +-----------+  +----------------+  |   |
|  |  | NV-Seg-CT |  | Nemotron  |  | NV-Gen-CT |  | NV-Gen-MR      |  |   |
|  |  | Client    |  | Nano Clt  |  | Client    |  | Client         |  |   |
|  |  +-----------+  +-----------+  +-----------+  +----------------+  |   |
|  |  +-----------+                                                     |   |
|  |  | NV-Reason |  (stub)                                            |   |
|  |  | CXR Client|                                                    |   |
|  |  +-----------+                                                     |   |
|  |                                                                    |   |
|  |  All 9 inherit BaseNIMClient: health check + retry + mock fallback|   |
|  +------------------------------------------------------------------+   |
|                                                                          |
|  +------------------------------------------------------------------+   |
|  | Workflow Engine (src/workflows/)                                   |   |
|  |                                                                    |   |
|  |  +-- CTHeadHemorrhageWorkflow    (< 90 sec, 3D U-Net)            |   |
|  |  +-- CTChestLungNoduleWorkflow   (< 5 min, RetinaNet+SegResNet)  |   |
|  |  +-- CTCoronaryAngiographyWorkflow (< 5 min, CAD-RADS)          |   |
|  |  +-- CXRRapidFindingsWorkflow    (< 30 sec, DenseNet-121)        |   |
|  |  +-- MRIBrainMSLesionWorkflow    (< 5 min, 3D U-Net+SyN)        |   |
|  |  +-- MRIProstateWorkflow         (< 5 min, PI-RADS v2.1)        |   |
|  |  +-- BreastBIRADSWorkflow        (< 5 min, BI-RADS)             |   |
|  |  +-- ThyroidTIRADSWorkflow       (< 3 min, TI-RADS)             |   |
|  |  +-- LiverLIRADSWorkflow         (< 5 min, LI-RADS)             |   |
|  |                                                                    |   |
|  |  All 9 inherit BaseImagingWorkflow: preprocess->infer->postprocess|   |
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
|  | 13 imaging collections     |   |                             |        |
|  | + 1 read-only genomic      |   | - Schema definitions        |        |
|  | 38,028 vectors total       |   | - CRUD operations           |        |
|  | IVF_FLAT / COSINE / 384d   |   |                             |        |
|  +----------------------------+   | - Parallel search            |        |
|                                   +-----------------------------+        |
|  +----------------------------+   +-----------------------------+        |
|  | Pydantic Models            |   | Ingest Pipelines            |        |
|  | src/models.py              |   | src/ingest/                 |        |
|  |                            |   |                             |        |
|  | - 13 collection models     |   | - PubMed parser             |        |
|  | - 9 NIM result models      |   | - ClinicalTrials parser     |        |
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
       |-- imaging_radiomics     (weight 0.08, top-5) --> 3 hits
       |-- imaging_reports       (weight 0.06, top-5) --> 4 hits
       |-- ... (all 14 collections including genomic_evidence)
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
    +-- LLMClient
    |       complete(messages) -> str
    |       (OpenAI-compatible /v1/chat/completions)
    |
    +-- NVSegmentCTClient
    |       segment_ct(input_path) -> SegmentationResult
    |
    +-- NemotronNanoClient
    |       reason(messages) -> str
    |
    +-- NVGenerateCTClient
    |       generate_ct(params) -> SyntheticCTResult
    |
    +-- NVGenerateMRClient
    |       generate_mr(params) -> SyntheticMRResult
    |
    +-- NVReasonCXRClient (stub)
            reason_cxr(image_path, question) -> ReasoningResult
```

### 5.3 NIMServiceManager

Coordinates all 9 NIM clients:

```python
NIMServiceManager(settings)
    .vista3d         -> VISTA3DClient
    .maisi           -> MAISIClient
    .vilam3          -> VILAM3Client
    .llm             -> LLMClient
    .nv_segment_ct   -> NVSegmentCTClient
    .nemotron_nano   -> NemotronNanoClient
    .nv_generate_ct  -> NVGenerateCTClient
    .nv_generate_mr  -> NVGenerateMRClient
    .nv_reason_cxr   -> NVReasonCXRClient (stub)
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

All nine reference workflows inherit from the same abstract base class:

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
    "ct_head_hemorrhage":     CTHeadHemorrhageWorkflow,
    "ct_chest_lung_nodule":   CTChestLungNoduleWorkflow,
    "ct_coronary_angiography": CTCoronaryAngiographyWorkflow,
    "cxr_rapid_findings":     CXRRapidFindingsWorkflow,
    "mri_brain_ms_lesion":    MRIBrainMSLesionWorkflow,
    "mri_prostate_pirads":    MRIProstateWorkflow,
    "breast_birads":          BreastBIRADSWorkflow,
    "thyroid_tirads":         ThyroidTIRADSWorkflow,
    "liver_lirads":           LiverLIRADSWorkflow,
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

**8 trigger conditions including:**
- Lung-RADS 4A+ findings --> queries EGFR, ALK, ROS1, KRAS variants
- CXR urgent consolidation --> queries infection-related genomic variants
- CAD-RADS >= 3 --> queries LDLR, PCSK9, APOB cardiovascular variants
- PI-RADS >= 4 --> queries BRCA2, HOXB13 cancer susceptibility
- BI-RADS 4+ --> queries BRCA1, BRCA2, ATM breast cancer variants
- TI-RADS TR4+ --> queries RET, BRAF thyroid cancer variants
- LI-RADS LR-4+ --> queries HFE, SERPINA1 liver disease variants
- Brain lesion high activity --> queries HLA-DRB1, MS susceptibility genes

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

### 10.2 Export Architecture (5 Formats)

The export system supports 5 formats: Markdown, JSON, PDF (ReportLab), FHIR R4, and DICOM SR (Structured Report via highdicom TID 1500). The `export_fhir()` function generates FHIR R4 DiagnosticReport Bundles with 103 SNOMED codes:

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

---

## 11. New Architectural Components (v2.1)

### 11.1 Agentic Reasoning (AIQ Toolkit)

The Clinical Imaging Engine integrates AIQ Plan/Execute/Reflect/Refine agentic reasoning with 6 tools for multi-step clinical analysis. The agent plans a series of tool invocations, executes them, reflects on intermediate results, and refines its approach before synthesizing a final answer.

### 11.2 NeMo Guardrails

NeMo Guardrails enforce PII protection (detecting and redacting patient identifiers), evidence grounding (ensuring claims are traceable to retrieved evidence), and disclaimer enforcement (appending clinical disclaimer to all outputs).

### 11.3 Radiomics (PyRadiomics-CUDA)

~1,500 radiomics features are extracted via PyRadiomics-CUDA, stored in the `imaging_radiomics` collection, and searchable via the RAG engine. Features include shape, first-order, GLCM, GLRLM, GLSZM, NGTDM, and GLDM descriptors.

### 11.4 Radiology Report NLP

A full radiology report parsing pipeline extracts findings, impressions, measurements, and coded diagnoses from free-text reports, storing structured results in the `imaging_reports` collection.

### 11.5 Protocol Optimization

12 ACR indications with patient-specific safety parameters. Protocol recommendations consider patient age, weight, renal function, contrast allergy history, and pregnancy status.

### 11.6 Dose Tracking

DRL (Diagnostic Reference Level) comparison with cumulative dose alerts. Tracks patient radiation exposure history and alerts when cumulative doses approach institutional thresholds.

### 11.7 Population Analytics (RAPIDS)

GPU-accelerated RAPIDS population analytics for cohort-level imaging trends, disease prevalence monitoring, and outcomes tracking across institutional imaging archives.

### 11.8 Streaming (Holoscan)

Holoscan real-time streaming pipeline for ultrasound and endoscopy, enabling sub-second AI inference on live video feeds.

### 11.9 MONAI Deploy MAPs

9 MONAI Application Packages (MAPs) packaged for clinical deployment, following MONAI Deploy standards for containerized inference pipelines.

### 11.10 MONAI Label

Interactive annotation with FLARE bridge, enabling radiologists to interactively segment structures and feed corrections back to the model training loop.

### 11.11 3D Visualization

Three.js rotating point cloud visualization for 3D volumetric data display in the React portal.

### 11.12 React Portal

**Live Analysis Layer.** A DICOMAnalyzer class (`src/dicom_analyzer.py`) provides real GPU inference on uploaded DICOM files. It auto-detects modality from DICOM headers, routes to the appropriate workflow, and runs actual model inference (DenseNet-121 for CXR, threshold segmentation for CT/MRI). Six MONAI model bundles (1.87 GB) are downloaded for production inference. API endpoints at `/analyze/*` handle file upload, sample analysis, and status reporting. The React portal exposes this at `/live-analysis` with drag-and-drop upload.

Full React portal with 10 pages, providing a modern web interface alongside the Streamlit workbench.

---

*For NIM-specific setup instructions, see `NIM_INTEGRATION_GUIDE.md`. For the complete implementation specification, see `PROJECT_BIBLE.md`.*
