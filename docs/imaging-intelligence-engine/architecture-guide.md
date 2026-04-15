# Clinical Imaging Engine -- Architecture Guide

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [System Diagram](#1-system-diagram)
2. [Component Interactions](#2-component-interactions)
3. [Data Flow](#3-data-flow)
4. [Collection Design Rationale](#4-collection-design-rationale)
5. [NIM Client Layer](#5-nim-client-layer)
6. [Workflow Pipeline](#6-workflow-pipeline)
7. [Cross-Modal Engine](#7-cross-modal-engine)
8. [Query Expansion](#8-query-expansion)
9. [RAG Pipeline](#9-rag-pipeline)
10. [Agent Orchestrator](#10-agent-orchestrator)
11. [Data Model Architecture](#11-data-model-architecture)

---

## 1. System Diagram

### 1.1 Full System Architecture

```
                          EXTERNAL USERS
                               |
                    +----------+----------+
                    |                     |
              +-----+------+      +------+-----+
              | Streamlit  |      | REST API   |
              | UI :8525   |      | :8524      |
              +-----+------+      +------+-----+
                    |                     |
                    +----------+----------+
                               |
                    +----------+----------+
                    |   Imaging Agent     |
                    |   (src/agent.py)    |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Query       |    | NIM Service   |    | Workflow      |
   | Expansion   |    | Manager       |    | Engine        |
   +------+------+    +-------+-------+    +-------+-------+
          |                    |                    |
          |           +-------+-------+            |
          |           |       |       |            |
          |        VISTA-3D MAISI  VILA-M3         |
          |        :8530  :8531  :8532             |
          |           |       |       |            |
          |           +-------+-------+            |
          |                    |                    |
          +--------------------+--------------------+
                               |
                    +----------+----------+
                    |    RAG Engine       |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Knowledge   |    | Milvus        |    | LLM           |
   | Graph       |    | Vector DB     |    | (Claude 4.6)  |
   | 25 pathols  |    | 11 Collections|    | Llama-3 NIM   |
   | 9 modalities|    |               |    | fallback      |
   | 21 anatomy  |    |               |    |               |
   +-------------+    +-------+-------+    +---------------+
                               |
                    +----------+----------+
                    |  etcd    |  MinIO   |
                    +----------+----------+
```

### 1.2 NVIDIA NIM Integration

```
  +------------------------------------------+
  |           NVIDIA NIM Microservices         |
  +------------------------------------------+
  | Vista3D   | MAISI     | VILA-M3  | Llama-3|
  | :8530     | :8531     | :8532    | :8520  |
  | 3D organ  | Synthetic | Visual   | Text   |
  | segment   | CT gen    | Q&A      | LLM    |
  +----+------+-----+-----+----+-----+---+---+
       |            |          |          |
       +------------+----------+----------+
                         |
                  +------v------+
                  | NIM Service |
                  | Manager     |
                  | (mock/live) |
                  +------+------+
                         |
                  +------v------+
                  | BaseNIMClient|
                  | health_check|
                  | retry logic |
                  | mock fallbk |
                  +-------------+
```

---

## 2. Component Interactions

### 2.1 Component Dependency Graph

```
ImagingUI (Streamlit) ──> FastAPI Server ──> ImagingAgent
                                               |
                                 +-------------+-------------+
                                 |             |             |
                          QueryExpansion  NIMServiceMgr   WorkflowEngine
                                 |             |             |
                                 +------+------+      6 Ref Workflows
                                        |                    |
                                   RAGEngine          CrossModalEngine
                                        |                    |
                                 +------+------+             |
                                 |             |             |
                            Milvus DB    Knowledge      genomic_evidence
                            (11 cols)     Graph          (shared col)
```

### 2.2 Module Responsibilities

| Module | File | Responsibilities |
|--------|------|-----------------|
| **Imaging Agent** | `src/agent.py` | Query classification, workflow dispatch, NIM orchestration, response assembly |
| **RAG Engine** | `src/rag_engine.py` | Multi-collection search, weighted scoring, query expansion, comparative analysis, LLM synthesis |
| **Knowledge Graph** | `src/knowledge.py` | 25 pathologies, 9 modalities, 21 anatomy entries |
| **Query Expansion** | `src/query_expansion.py` | Domain-specific maps, keyword-to-term expansion, entity resolution |
| **NIM Service Manager** | `src/nim/service_manager.py` | VISTA-3D, MAISI, VILA-M3, Llama-3 client coordination |
| **Workflow Engine** | `src/workflows/` | 6 reference workflows (CT head, CT chest, CT coronary, CXR, MRI brain, MRI prostate) |
| **Collections** | `src/collections.py` | Schema definitions, CRUD operations, parallel search |
| **Cross-Modal** | `src/cross_modal.py` | Imaging-to-genomics trigger evaluation |
| **Export** | `src/export.py` | Markdown, JSON, PDF, FHIR R4 DiagnosticReport |
| **Models** | `src/models.py` | 10 collection models, 4 NIM result models, search models, agent I/O |
| **Ingest** | `src/ingest/` | PubMed parser, ClinicalTrials parser, 6 seed data parsers |

### 2.3 Interface Contracts

**Imaging Agent inputs/outputs:**
```
Input:  ImagingQuery(question, modality_filter?, body_region?, patient_context?)
Output: ImagingResponse(answer, evidence, nim_status, workflow_results, cross_modal)
```

**Workflow Engine inputs/outputs:**
```
Input:  WorkflowRequest(workflow_name, input_path?, study_data?)
Output: WorkflowResult(findings, measurements, classification, severity, inference_time_ms, is_mock)
```

**NIM Service Manager inputs/outputs:**
```
Input:  NIMRequest(endpoint, payload, client_type)
Output: NIMResponse(result, status, latency_ms, is_mock)
```

---

## 3. Data Flow

### 3.1 RAG Query Flow

```
Step 1: RECEIVE QUERY
  User Query: "What is ACR Lung-RADS classification?"
  |
Step 2: QUERY CLASSIFICATION
  Detect comparative ("X vs Y")? --> No
  Detect modality filter? --> CT
  Detect body region? --> Chest
  |
Step 3: QUERY EXPANSION
  "Lung-RADS" --> ["lung_rads", "lung_cancer_screening",
                    "nodule_management", "ACR", ...]
  |
Step 4: EMBEDDING
  BGE-small-en-v1.5: "Represent this sentence: ..."
  Output: 384-dim float32 vector
  |
Step 5: PARALLEL MULTI-COLLECTION SEARCH
  imaging_literature    (weight 0.18, top-5) --> 5 hits
  imaging_findings      (weight 0.15, top-5) --> 3 hits
  imaging_guidelines    (weight 0.10, top-5) --> 5 hits
  imaging_trials        (weight 0.12, top-5) --> 4 hits
  ... (all 11 collections)
  |
Step 6: WEIGHTED SCORE MERGE
  Combine hits across collections
  Apply collection weights
  Filter by SCORE_THRESHOLD (0.4)
  Sort by weighted score descending
  |
Step 7: KNOWLEDGE GRAPH AUGMENTATION
  Match "lung_nodule" pathology entry
  Inject: Lung-RADS categories, severity criteria,
          CT characteristics, AI models
  |
Step 8: LLM SYNTHESIS
  Build prompt: question + evidence + knowledge context
  Inject conversation history (up to 3 prior turns)
  Call Claude API (or Llama-3 NIM fallback)
  |
Step 9: RESPONSE ASSEMBLY
  Grounded answer with evidence citations
  Source references with scores
  Follow-up question suggestions
  NIM service availability status
```

### 3.2 Workflow Execution Flow

```
API Request: POST /workflow/ct_head_hemorrhage/run
  |
  v
[1. Workflow Registry Lookup]
  WORKFLOW_REGISTRY["ct_head_hemorrhage"]
  Instantiate CTHeadHemorrhageWorkflow(mock_mode=True)
  |
  v
[2. Preprocess]
  Mock: skip (return synthetic volume metadata)
  Live: LoadImaged -> EnsureChannelFirst -> Orientationd(RAS)
        -> Spacingd(1mm) -> ScaleIntensityRanged(0-80 HU)
  |
  v
[3. Infer]
  Mock: return synthetic segmentation result
  Live: 3D U-Net binary segmentation via MONAI
  |
  v
[4. Postprocess]
  Volume estimation: voxel count x voxel volume
  Midline shift: center of mass vs falx cerebri
  Max thickness measurement
  BTF urgency classification (P1/P2/P4)
  |
  v
[5. WorkflowResult]
  findings: [{category, description, severity, recommendation}]
  measurements: {volume_ml, shift_mm, thickness_mm}
  classification: "P1" / "P2" / "P4"
  severity: critical / urgent / routine
  inference_time_ms, is_mock
```

### 3.3 DICOM Ingestion Flow

```
Orthanc DICOM Server (port 8042 HTTP, 4242 C-STORE)
    |
    v
POST /events/dicom-webhook (study.complete event)
    |
    v
determine_workflow(modality, body_region) -> workflow name
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

---

## 4. Collection Design Rationale

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

### 4.2 Collection Inventory

| # | Collection Name | Description | Weight | Seed Records |
|---|----------------|-------------|--------|-------------|
| 1 | imaging_literature | PubMed research papers, radiology reviews | 0.18 | 464 |
| 2 | imaging_findings | Imaging finding templates and clinical patterns | 0.15 | 80 |
| 3 | imaging_trials | Clinical trials for imaging interventions | 0.12 | 60 |
| 4 | imaging_guidelines | ACR, RSNA, ESR guideline recommendations | 0.10 | 50 |
| 5 | imaging_protocols | Imaging protocol parameters by modality | 0.09 | 45 |
| 6 | imaging_anatomy | Anatomical reference data and landmarks | 0.08 | 40 |
| 7 | imaging_devices | FDA-cleared AI devices and 510(k) database | 0.07 | 35 |
| 8 | imaging_benchmarks | AI model benchmark results and metrics | 0.06 | 30 |
| 9 | imaging_report_templates | Structured reporting templates (BI-RADS, Lung-RADS, etc.) | 0.06 | 40 |
| 10 | imaging_pathology_correlation | Imaging-pathology correlation data | 0.05 | 32 |
| 11 | genomic_evidence | Shared VCF-derived genomic variants (read-only) | 0.04 | -- |
| | **Total** | | **1.00** | **876** |

### 4.3 Collection Weight Distribution

```
imaging_literature      ██████████████████████  0.18
imaging_findings        ██████████████████      0.15
imaging_trials          ████████████████        0.12
imaging_guidelines      █████████████           0.10
imaging_protocols       ████████████            0.09
imaging_anatomy         ██████████              0.08
imaging_devices         █████████               0.07
imaging_benchmarks      ████████                0.06
imaging_report_templ.   ████████                0.06
imaging_path_corr.      ███████                 0.05
genomic_evidence        █████                   0.04
                                           Sum: 1.00
```

### 4.4 Schema Pattern

Every collection follows the same field pattern:

```python
FieldSchema("id",        DataType.VARCHAR, max_length=100, is_primary=True)
FieldSchema("embedding", DataType.FLOAT_VECTOR, dim=384)
FieldSchema("text",      DataType.VARCHAR, max_length=3000)
# ... domain-specific metadata fields (VARCHAR, INT64, FLOAT, etc.)
```

### 4.5 Search Strategy

1. **Parallel search:** All collections searched simultaneously using ThreadPoolExecutor
2. **Per-collection top-K:** Default 5 results per collection (configurable 1-50)
3. **Weighted scoring:** Each collection has a configurable weight (0.04 to 0.18)
4. **Score threshold:** Results below 0.4 cosine similarity are filtered out
5. **Asymmetric embedding:** Queries use BGE instruction prefix

---

## 5. NIM Client Layer

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
    |       Port: 8530
    |
    +-- MAISIClient
    |       generate(body_region, resolution) -> SyntheticCTResult
    |       Port: 8531
    |
    +-- VILAM3Client
    |       analyze(question, input_path) -> VLMResponse
    |       Port: 8532
    |
    +-- LlamaLLMClient
            complete(messages) -> str
            Port: 8520 (OpenAI-compatible /v1/chat/completions)
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

Every NIM-dependent feature degrades gracefully through three tiers: Full (local GPU) -> Cloud (API endpoint) -> Mock (synthetic results). This ensures the agent remains functional for RAG queries and knowledge retrieval even without GPU hardware.

---

## 6. Workflow Pipeline

### 6.1 BaseImagingWorkflow (ABC)

All six reference workflows inherit from the same abstract base class:

```python
class BaseImagingWorkflow(ABC):
    WORKFLOW_NAME: str
    TARGET_LATENCY_SEC: float
    MODALITY: str
    BODY_REGION: str
    MODELS_USED: List[str]

    preprocess(input_path)    -> Any
    infer(preprocessed)       -> Dict
    postprocess(result)       -> WorkflowResult
    _mock_inference()         -> Dict
    run(input_path)           -> WorkflowResult
    get_workflow_info()       -> Dict
```

### 6.2 Workflow Registry

```python
WORKFLOW_REGISTRY = {
    "ct_head_hemorrhage":    CTHeadHemorrhageWorkflow,     # < 90 sec, 3D U-Net
    "ct_chest_lung_nodule":  CTChestLungNoduleWorkflow,    # < 5 min, RetinaNet+SegResNet
    "ct_coronary_angiography": CTCoronaryAngiographyWorkflow,  # Calcium scoring, stenosis
    "cxr_rapid_findings":    CXRRapidFindingsWorkflow,     # < 30 sec, DenseNet-121
    "mri_brain_ms_lesion":   MRIBrainMSLesionWorkflow,     # < 5 min, 3D U-Net+SyN
    "mri_prostate_pirads":   MRIProstateWorkflow,          # PI-RADS classification
}
```

Dynamic dispatch via the `/workflow/{name}/run` API endpoint.

### 6.3 Workflow Details

| Workflow | Modality | Region | Model | Output |
|----------|----------|--------|-------|--------|
| CT Head Hemorrhage | CT | Head | SegResNet 3D U-Net | Volume (mL), midline shift (mm), urgency (P1/P2/P4) |
| CT Chest Lung Nodule | CT | Chest | RetinaNet + SegResNet | Nodule count, volume, doubling time, Lung-RADS |
| CT Coronary Angiography | CT | Chest | Calcium scoring | Agatston score, stenosis grade, CAD-RADS |
| CXR Rapid Findings | CXR | Chest | DenseNet-121 (CheXpert) | Multi-label classification, confidence scores |
| MRI Brain MS Lesion | MRI | Brain | UNEST + SyN registration | Lesion count, volume, longitudinal change |
| MRI Prostate PI-RADS | MRI | Pelvis | Lesion detection | PI-RADS score, lesion locations, volume |

### 6.4 Error Handling

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

## 7. Cross-Modal Engine

### 7.1 Genomic Pipeline Trigger

The CrossModalTrigger class (`src/cross_modal.py`) automatically enriches high-risk imaging findings with genomic context from the shared `genomic_evidence` collection (3.5M vectors from Stage 2).

**Trigger conditions:**
- Lung-RADS 4A+ findings --> queries EGFR, ALK, ROS1, KRAS variants
- CXR urgent consolidation --> queries infection-related genomic variants
- Brain lesion with mass effect --> queries tumor-associated variants

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

### 7.2 FHIR R4 Export Architecture

The `export_fhir()` function generates FHIR R4 DiagnosticReport Bundles:

```
FHIR Bundle (type: collection)
+-- Patient resource (stub with identifier)
+-- ImagingStudy resource (modality auto-detected from query)
+-- Observation resources (one per workflow finding)
|   +-- SNOMED CT coding (finding category)
|   +-- Interpretation (severity -> HH/H/A/N)
|   +-- Components (measurements with UCUM units)
+-- DiagnosticReport resource
    +-- LOINC category (LP29684-5 Radiology)
    +-- LOINC code (18748-4 Diagnostic imaging study)
    +-- conclusionCode (SNOMED for all findings)
    +-- extension (cross-modal enrichment summary)
```

### 7.3 Drug Discovery Pipeline Feed

```
Quantitative imaging endpoint
    |-- Tumor volume change
    |-- RECIST measurements
    |-- Treatment response
    |
    v
Drug Discovery Pipeline (Stage 3)
    |-- Treatment-response tracking
    |-- Molecular target validation
```

---

## 8. Query Expansion

### 8.1 Domain-Specific Maps

The query expansion module resolves imaging-specific terminology:

| Category | Example Input | Example Expansions |
|----------|---------------|-------------------|
| Modality | CT | computed tomography, MDCT, helical CT |
| Body region | chest | thorax, thoracic, pulmonary |
| Finding | nodule | lung nodule, SPN, pulmonary nodule, GGO |
| Classification | Lung-RADS | lung_rads, lung_cancer_screening, ACR |
| AI model | VISTA-3D | vista3d, 3D segmentation, organ segmentation |
| Measurement | SUV | standardized uptake value, FDG uptake |
| Protocol | contrast | contrast-enhanced, IV contrast, iodinated |
| Anatomy | coronary | coronary artery, LAD, RCA, circumflex |

### 8.2 Entity Resolution

Abbreviations and synonyms are resolved before embedding:

- "CXR" -> "chest X-ray, chest radiograph"
- "MRI" -> "magnetic resonance imaging"
- "GGO" -> "ground-glass opacity"
- "LGE" -> "late gadolinium enhancement"
- "PE" -> "pulmonary embolism"

---

## 9. RAG Pipeline

### 9.1 RAG Architecture

```
SearchPlan (from Query Expansion)
        |
  Embedding (BGE-small-en-v1.5, 384-dim)
        |
  Multi-Collection Search (11 collections, parallel)
        |
  Score Filtering (threshold 0.4)
        |
  Weight Application (per-collection weights)
        |
  Deduplication (content hash)
        |
  Citation Scoring
  - High confidence: score >= 0.80
  - Medium confidence: score >= 0.60
  - Standard: score < 0.60
        |
  Knowledge Graph Augmentation
  - 25 pathology entries
  - 9 modality entries
  - 21 anatomy entries
        |
  LLM Synthesis
  - Primary: Claude Sonnet 4.6
  - Fallback: Llama-3 8B NIM (port 8520)
  - System prompt: Radiology domain expert
  - Conversation memory: 3 turns
        |
  Response Assembly
  - Grounded answer with citations
  - NIM service status
  - Follow-up suggestions
```

### 9.2 Embedding Configuration

| Parameter | Value |
|-----------|-------|
| Model | BAAI/bge-small-en-v1.5 |
| Parameters | 33M |
| Dimensions | 384 |
| Metric | COSINE |
| Index type | IVF_FLAT (nlist=1024, nprobe=16) |
| Batch size | 32 |
| Runtime | CPU (no GPU required) |
| Instruction prefix | "Represent this sentence for searching relevant passages: " |
| Search mode | Asymmetric |

---

## 10. Agent Orchestrator

### 10.1 Key Statistics

| Metric | Value |
|--------|-------|
| Reference workflows | 6 (CT head, CT chest, CT coronary, CXR, MRI brain, MRI prostate) |
| Demo cases | 4 (DEMO-001 through DEMO-004) |
| Milvus collections | 11 (10 imaging + 1 read-only genomic_evidence) |
| Seed imaging vectors | 876 across 10 owned collections |
| NVIDIA NIMs | 3 (Vista3D, MAISI, VILA-M3) + Llama-3 fallback |
| Streamlit tabs | 9 |
| Docker Compose services | 13 full mode, 6 lite mode |
| Unit tests | 539, E2E checks 9/9 |
| Output formats | 4 (Markdown, JSON, PDF, FHIR R4) |

### 10.2 Error Handling Strategy

The agent implements graceful degradation at every layer:

1. **NIM unavailable**: Falls back to mock mode automatically (Full -> Cloud -> Mock)
2. **Milvus unavailable**: Returns error with clear message
3. **LLM unavailable (Claude)**: Falls back to Llama-3 NIM, then to mock
4. **Workflow error**: Returns WorkflowResult with status=FAILED
5. **Cross-modal error**: Logs warning, returns response without genomic enrichment
6. **Timeout**: Returns partial results with timeout indicator

### 10.3 API Endpoints Summary

| Router | Prefix | Endpoints |
|--------|--------|-----------|
| meta_agent | `/api` | `/api/ask` |
| nim | `/nim` | `/nim/status`, `/nim/vista3d/segment`, `/nim/maisi/generate`, `/nim/vilam3/analyze` |
| workflows | (root) | `/workflows`, `/workflow/{name}/info`, `/workflow/{name}/run` |
| reports | (root) | `/reports/generate` |
| events | `/events` | `/events/dicom-webhook`, `/events/history`, `/events/status` |
| core | (root) | `/health`, `/collections`, `/query`, `/search`, `/find-related`, `/knowledge/stats`, `/metrics` |

### 10.4 Port Configuration

| Service | Port | Protocol |
|---------|------|----------|
| FastAPI REST API | 8524 | HTTP |
| Streamlit UI | 8525 | HTTP |
| VISTA-3D NIM | 8530 | HTTP |
| MAISI NIM | 8531 | HTTP |
| VILA-M3 NIM | 8532 | HTTP |
| Llama-3 NIM | 8520 | HTTP |
| Milvus (shared) | 19530 | gRPC |

---

## 11. Data Model Architecture

### 11.1 Model Hierarchy

```
                    ImagingResponse (top-level output)
                    /        |         \           \
              answer    evidence    nim_status    workflow_results
                |           |            |              |
              str    List[Dict]    Dict[str,str]  List[WorkflowResult]

                    WorkflowResult (per-workflow output)
                    /        |         \           \
              findings  measurements  classification  severity
                |           |            |              |
           List[Finding]  Dict      str           SeverityLevel
```

### 11.2 Ingest Pipeline

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

### 11.3 Scheduling

APScheduler (`src/scheduler.py`) supports periodic re-ingestion:
- Default interval: 168 hours (weekly)
- Configurable via `IMAGING_INGEST_SCHEDULE_HOURS`
- Disabled by default (`IMAGING_INGEST_ENABLED=false`)

---

*For NIM-specific setup instructions, see `NIM_INTEGRATION_GUIDE.md`. For the complete implementation specification, see `PROJECT_BIBLE.md`.*

---

!!! warning "Clinical Decision Support Disclaimer"
    The Clinical Imaging Engine is a clinical decision support research tool for medical image analysis. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
