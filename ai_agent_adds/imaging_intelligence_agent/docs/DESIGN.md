# Imaging Intelligence Agent -- Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 2.0.0
**License:** Apache 2.0

---

## 1. Executive Summary

### 1.1 Platform Context: HCLS AI Factory 3-Engine Architecture

The Imaging Intelligence Agent operates within the HCLS AI Factory, a three-engine precision medicine platform on NVIDIA DGX Spark:

- **Stage 1 -- Genomics Engine:** Parabricks/DeepVariant/BWA-MEM2 for FASTQ-to-VCF variant calling (GPU-accelerated)
- **Stage 2 -- RAG/Chat Engine:** Milvus (3.56M vectors) + Claude AI for variant interpretation and evidence synthesis
- **Stage 3 -- Drug Discovery Engine:** BioNeMo MolMIM/DiffDock/RDKit for lead optimization across 171 druggable targets

The platform includes **11 intelligence agents**:

| # | Agent | Port | Domain |
|---|---|---|---|
| 1 | Biomarker Intelligence | :8529 | Biomarker discovery and stratification |
| 2 | Oncology Intelligence | :8527/:8528 | Cancer genomics and treatment selection |
| 3 | CAR-T Intelligence | -- | CAR-T cell therapy development |
| 4 | **Imaging Intelligence** | **:8524** | **Medical imaging AI (this agent)** |
| 5 | Autoimmune Intelligence | -- | Autoimmune disease genomics |
| 6 | Pharmacogenomics Intelligence | :8107 | Drug metabolism and dosing |
| 7 | Clinical Trial Intelligence | :8538 | Trial design and patient matching |
| 8 | Rare Disease Diagnostic | :8134 | Rare disease diagnosis |
| 9 | Single-Cell Intelligence | :8540 | Single-cell transcriptomics |
| 10 | Cardiology Intelligence | :8126 | Cardiac genetics and risk |
| 11 | Neurology Intelligence | -- | Neurological genetics |

### 1.2 Agent Summary

The Imaging Intelligence Agent extends the HCLS AI Factory platform to support automated detection, segmentation, longitudinal tracking, and clinical triage of medical imaging studies. The agent processes CT, MRI, and chest X-ray studies on DGX Spark hardware using a multi-collection RAG engine backed by Milvus 2.4, four NVIDIA NIM microservices, and six reference clinical workflows -- from DICOM ingestion through agentic inference to structured clinical output.

Six reference workflows cover the highest-impact radiology use cases:

1. **CT Head Hemorrhage Triage** -- SegResNet segmentation, volume estimation, midline shift, urgency routing
2. **CT Chest Lung Nodule Tracking** -- RetinaNet detection + SegResNet segmentation, volume doubling time, Lung-RADS
3. **CT Coronary Angiography** -- Calcium scoring, stenosis grading, CAD-RADS classification
4. **CXR Rapid Findings** -- DenseNet-121 multi-label classification (CheXpert pretrained weights)
5. **MRI Brain MS Lesion Tracking** -- UNEST segmentation, longitudinal lesion matching
6. **MRI Prostate PI-RADS** -- PI-RADS scoring and lesion characterization

Cross-modal triggers link imaging to genomics: Lung-RADS 4A+ findings automatically query the shared genomic_evidence collection (3.5M vectors from Stage 2) for relevant genetic variants (EGFR, ALK, ROS1, KRAS).

### Key Metrics

| Metric | Value |
|---|---|
| Reference workflows | **6** (CT head, CT chest, CT cardiac, CXR, MRI brain, MRI prostate) |
| Milvus collections | **11** (10 imaging + 1 read-only genomic_evidence) |
| Seed vectors | **876** across 10 imaging collections |
| Embedding model | **BAAI/bge-small-en-v1.5** (384-dim, IVF_FLAT, COSINE) |
| LLM | **Claude Sonnet 4.6** (Anthropic) primary, **Llama-3 8B NIM** fallback |
| Docker Compose services | **13** full mode, **6** lite mode |
| Output formats | **4** (Markdown, JSON, PDF, FHIR R4 DiagnosticReport) |
| Demo cases | **4** pre-loaded (DEMO-001 through DEMO-004) |
| UI tabs | **9** (Evidence Explorer, Workflow Runner, Image Gallery, Protocol Advisor, Device & AI Ecosystem, Dose Intelligence, Reports & Export, Patient 360, Benchmarks & Validation) |

---

## 2. Architecture Overview

### 2.1 Design Principles

1. **Graceful degradation** -- Every NIM-dependent feature falls back automatically (Full -> Cloud -> Mock)
2. **Shared infrastructure** -- Reuses Milvus 2.4 and BGE-small-en-v1.5 from the HCLS AI Factory platform
3. **Cross-agent interoperability** -- Read-only access to `genomic_evidence` collection from Stage 2
4. **Consistent patterns** -- Follows the same Pydantic BaseSettings, collection manager, and RAG engine patterns as the other intelligence agents

### 2.2 Component Diagram

```
+=========================================================================+
|  PRESENTATION:  Streamlit Imaging Workbench (8525)  |  FastAPI (8524)   |
|                 9-tab UI with guided tour            |  5 route modules  |
+=========================================================================+
                    |                            |
+=========================================================================+
|  INTELLIGENCE:   Imaging Agent         |  RAG Engine                    |
|                  src/agent.py          |  src/rag_engine.py             |
|                  Query classification  |  Multi-collection search       |
|                  Workflow dispatch      |  Weighted scoring + expansion  |
|                  Response assembly     |  LLM synthesis (Claude/Llama)  |
|                                        |                                |
|                  Knowledge Graph       |  Query Expansion               |
|                  src/knowledge.py      |  src/query_expansion.py        |
|                  25 pathologies         |  Domain-specific term maps     |
|                  9 modalities           |                                |
|                  21 anatomy regions     |                                |
+=========================================================================+
                    |                            |
+=========================================================================+
|  INFERENCE:      NIM Service Manager (src/nim/service_manager.py)       |
|                  VISTA-3D (8530) | MAISI (8531) | VILA-M3 (8532)       |
|                  Llama-3 8B (8520) -- all via BaseNIMClient ABC         |
|                                                                          |
|                  Workflow Engine (src/workflows/)                        |
|                  CTHead | CTChest | CTCoronary | CXR | MRIBrain |       |
|                  MRIProstate -- all via BaseImagingWorkflow              |
+=========================================================================+
                    |                            |
+=========================================================================+
|  DATA:           Milvus 2.4 (19530)   |  Ingest Pipelines              |
|                  10 imaging collections |  PubMed, ClinicalTrials.gov   |
|                  + 1 genomic (read-only)|  8 seed data parsers           |
|                  IVF_FLAT/COSINE/384d  |  APScheduler integration       |
+=========================================================================+
```

### 2.3 VAST AI OS Mapping

The Imaging Intelligence Agent maps to VAST AI OS infrastructure as follows:

| HCLS Component | VAST AI OS Resource | Notes |
|---|---|---|
| Milvus 2.4 vector store | VAST DataStore (NFS/S3) backing Milvus persistent volumes | `etcd_data`, `minio_data`, `milvus_data` volumes |
| BGE-small-en-v1.5 embeddings | VAST Compute (CPU) | Embedding inference on CPU, ~384-dim vectors |
| Claude Sonnet 4.6 LLM | External API (Anthropic) | Primary LLM, API key in `.env` |
| Llama-3 8B NIM | VAST Compute (GPU) | Fallback LLM, OpenAI-compatible endpoint |
| VISTA-3D / MAISI / VILA-M3 | VAST Compute (GPU) | NIM microservices, mock fallback when GPU unavailable |
| FastAPI + Streamlit | VAST Compute (CPU) | Application tier, ports 8524/8525 |
| Orthanc DICOM server | VAST DataStore (NFS) backing `orthanc_data` volume | DICOM C-STORE on 4242, REST on 8042 |
| Docker Compose orchestration | VAST Container Runtime | 13-service full stack or 6-service lite |
| Prometheus metrics | VAST Monitoring | `/metrics` endpoint on FastAPI |

VAST AI OS provides the unified storage layer for Milvus vector persistence, DICOM image storage, and NIM model weight caching. The GPU compute resources host the NIM inference containers while CPU resources handle the application tier.

---

## 3. Data Flow

### 3.1 RAG Query Flow

```
User Query --> Query Classification (modality, body region, comparative)
          --> Query Expansion (domain-specific term maps)
          --> BGE-small-en-v1.5 embedding (384-dim)
          --> Parallel multi-collection search (11 collections, ThreadPoolExecutor)
          --> Weighted score merge (per-collection weights 0.04-0.18, threshold 0.4)
          --> Knowledge graph augmentation (pathology context injection)
          --> LLM synthesis (Claude Sonnet 4.6 primary, Llama-3 NIM fallback)
          --> Response assembly (answer + citations + follow-ups + NIM status)
```

### 3.2 Workflow Execution Flow

```
POST /workflow/{name}/run
    --> WORKFLOW_REGISTRY lookup --> Instantiate workflow
    --> preprocess (MONAI transforms or mock skip)
    --> infer (model inference or mock synthetic result)
    --> postprocess (measurements, classification, severity)
    --> WorkflowResult (findings, measurements, severity, inference_time_ms, is_mock)
```

### 3.3 Demo Case Flow

```
POST /demo-cases/{case_id}/run
    --> Load pre-configured clinical scenario
    --> Execute appropriate workflow with mock data
    --> Enrich with cross-modal genomic context
    --> Return results with talking points
```

---

## 4. Milvus Collection Design

### 4.1 Index Configuration

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist / nprobe | 1024 / 16 |
| Dimension | 384 |
| Embedding model | BAAI/bge-small-en-v1.5 |

### 4.2 Collections

| Collection | Seed Records | Content |
|---|---|---|
| `imaging_literature` | 50 | PubMed research papers and reviews |
| `imaging_trials` | 40 | ClinicalTrials.gov AI-in-imaging records |
| `imaging_findings` | 50 | Imaging finding templates and patterns |
| `imaging_protocols` | 40 | Acquisition protocols and parameters |
| `imaging_devices` | 50 | FDA-cleared AI/ML medical devices |
| `imaging_anatomy` | 30 | Anatomical structure references |
| `imaging_benchmarks` | 40 | Model performance benchmarks |
| `imaging_guidelines` | 40 | Clinical practice guidelines (ACR, RSNA, NCCN) |
| `imaging_report_templates` | 50 | Structured radiology report templates |
| `imaging_datasets` | 50 | Public imaging datasets (TCIA, PhysioNet) |
| `genomic_evidence` | 3,561,170 | *(read-only)* Shared from Stage 2 RAG pipeline |

**Total: 876 seed vectors across 10 imaging collections + 3,561,170 genomic evidence vectors.**

### 4.3 Search Strategy

1. **Parallel search:** All collections searched simultaneously via `ThreadPoolExecutor`
2. **Per-collection top-K:** Default 5 results per collection (configurable 1-50)
3. **Weighted scoring:** Each collection has a configurable weight (0.04 to 0.18)
4. **Score threshold:** Results below 0.4 cosine similarity filtered out
5. **Asymmetric embedding:** Queries use BGE instruction prefix for retrieval

---

## 5. NIM Client Layer

### 5.1 Client Hierarchy

All NIM clients inherit from `BaseNIMClient` (ABC), which provides cached health checks (30s interval), exponential-backoff retry via tenacity (3 attempts), and automatic mock fallback.

```
BaseNIMClient (ABC)
    +-- VISTA3DClient       segment()         port 8530   132 anatomical classes
    +-- MAISIClient         generate()        port 8531   Synthetic CT volumes
    +-- VILAM3Client        analyze()         port 8532   Vision-language for radiology
    +-- LlamaLLMClient      complete()        port 8520   OpenAI-compatible chat
```

### 5.2 Fallback Logic

```
Request --> NIM available? --> Yes --> Real NIM response
                           --> No  --> Cloud NIM? --> Yes --> Cloud API response
                                                  --> No  --> Mock enabled? --> Yes --> Synthetic mock response
                                                                            --> No  --> ConnectionError
```

The `NIMServiceManager` coordinates all four clients and exposes `check_all_services()` returning status per NIM (available / mock / unavailable).

---

## 6. Workflow Pipeline

### 6.1 Reference Workflows

| Workflow | Modality | Latency | Model |
|---|---|---|---|
| CT Head Hemorrhage Triage | CT | < 90 sec | SegResNet (MONAI wholeBody_ct_segmentation) |
| CT Chest Lung Nodule Tracking | CT | < 5 min | RetinaNet + SegResNet (MONAI lung_nodule_ct_detection) |
| CT Coronary Angiography | CT | < 3 min | Calcium scoring + stenosis grading |
| CXR Rapid Findings | X-ray | < 30 sec | DenseNet-121 (torchxrayvision CheXpert) |
| MRI Brain MS Lesion Tracking | MRI | < 5 min | UNEST (MONAI wholeBrainSeg_Large_UNEST_segmentation) |
| MRI Prostate PI-RADS | MRI | < 5 min | PI-RADS scoring and lesion characterization |

### 6.2 BaseImagingWorkflow (ABC)

All workflows implement `preprocess() -> infer() -> postprocess() -> WorkflowResult` with full mock mode support. The `WORKFLOW_REGISTRY` dict enables dynamic dispatch via API:

```python
WORKFLOW_REGISTRY = {
    "ct_head_hemorrhage":     CTHeadHemorrhageWorkflow,
    "ct_chest_lung_nodule":   CTChestLungNoduleWorkflow,
    "ct_coronary_angiography": CTCoronaryAngiographyWorkflow,
    "cxr_rapid_findings":     CXRRapidFindingsWorkflow,
    "mri_brain_ms_lesion":    MRIBrainMSLesionWorkflow,
    "mri_prostate_pirads":    MRIProstatePIRADSWorkflow,
}
```

---

## 7. Ingest Pipelines

### 7.1 Pattern

```
[Source] --> fetch() --> parse() --> embed() --> store()
         HTTP/API    Extract     BGE-small   Milvus upsert
         PubMed      Normalize   384-dim
         CT.gov      Pydantic
         Seed JSON   Validate
```

### 7.2 Parsers

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
| `dicom_watcher.py` | Orthanc DICOM polling | DICOM event routing |

APScheduler (`src/scheduler.py`) supports weekly re-ingestion, disabled by default.

---

## 8. API and UI Layers

### 8.1 FastAPI (port 8524)

Lifespan-managed initialization of Milvus, embedding model, NIM services, and RAG engine. CORS enabled. Prometheus metrics at `/metrics`.

| Router | Prefix | Key Endpoints |
|---|---|---|
| `meta_agent` | `/api` | `/api/ask` |
| `nim` | `/nim` | `/nim/status`, `/nim/vista3d/segment`, `/nim/maisi/generate`, `/nim/vilam3/analyze` |
| `workflows` | (root) | `/workflows`, `/workflow/{name}/info`, `/workflow/{name}/run` |
| `reports` | (root) | `/reports/generate` |
| `events` | `/events` | `/events/dicom-webhook`, `/events/history`, `/events/status` |

Core endpoints: `/health`, `/collections`, `/query`, `/search`, `/find-related`, `/knowledge/stats`, `/metrics`, `/demo-cases`, `/demo-cases/{id}/run`, `/protocol/recommend`, `/dose/reference`, `/dose/summary`

### 8.2 Streamlit Imaging Workbench (port 8525)

9-tab interface with NVIDIA dark/green styling:

| Tab | Features |
|---|---|
| **Evidence Explorer** | Multi-collection RAG queries, pre-filled example queries, comparative analysis auto-detection |
| **Workflow Runner** | 6-step pipeline animation, annotated AI images, 4 demo case dropdown, mock/live mode |
| **Image Gallery** | 5 CXR pathologies showcase, cross-modality gallery, 3D volume slice viewer |
| **Protocol Advisor** | Patient-specific imaging protocol recommendations, pre-filled examples |
| **Device & AI Ecosystem** | 50 FDA-cleared/research AI devices, filterable by modality and task |
| **Dose Intelligence** | 20 protocols with AI dose reduction, Plotly charts, 36% avg reduction |
| **Reports & Export** | Markdown, JSON, NVIDIA-themed PDF, FHIR R4 DiagnosticReport export |
| **Patient 360** | Cross-modal genomic enrichment, interactive Plotly network graph |
| **Benchmarks & Validation** | Plotly donut chart, model performance benchmarks, collection stats |

Additional features: sidebar guided tour (9-step demo flow), NIM status indicators, multi-turn conversation memory.

---

## 9. Cross-Modal Integration

### 9.1 Genomic Pipeline Trigger

The `CrossModalTrigger` (`src/cross_modal.py`) enriches high-risk imaging findings with genomic context from the shared `genomic_evidence` collection (3.5M vectors).

- **Lung-RADS 4A+** --> queries EGFR, ALK, ROS1, KRAS variants
- **CXR urgent consolidation** --> queries infection-related genomic variants
- **Brain lesion high activity** --> queries HLA-DRB1, demyelination markers

### 9.2 FHIR R4 Export

`export_fhir()` generates FHIR R4 DiagnosticReport Bundles with Patient, ImagingStudy, Observation (per finding with SNOMED CT coding, UCUM measurements), and DiagnosticReport (LOINC category, cross-modal enrichment extension).

### 9.3 DICOM Ingestion

Orthanc webhook (POST `/events/dicom-webhook`) receives study.complete events, determines the appropriate workflow from modality and body region, executes it, and stores results in event history (in-memory, max 200 entries).

### 9.4 Export Formats

| Format | Function | Output |
|---|---|---|
| Markdown | `export_markdown()` | Structured clinical report |
| JSON | `export_json()` | Full Pydantic model serialization |
| PDF | `export_pdf()` | NVIDIA-themed report (ReportLab) |
| FHIR R4 | `export_fhir()` | DiagnosticReport Bundle (SNOMED CT, LOINC, DICOM) |

### 9.5 Cross-Agent Integration (4 Agents)

Beyond genomic cross-modal triggers, the Imaging Intelligence Agent integrates with 4 peer agents:

| Peer Agent | Port | Trigger | Integration Purpose |
|---|---|---|---|
| Oncology Agent | :8527 | Lung-RADS 4A+ or suspicious mass | Molecular profiling trigger, tumor staging correlation |
| Clinical Trial Agent | :8538 | Imaging-based trial eligibility | RECIST/RANO response criteria, imaging endpoint matching |
| Cardiology Agent | :8126 | CAD-RADS findings, cardiac CT | Coronary risk scoring, FH gene queries (LDLR, PCSK9, APOB) |
| Neurology Agent | -- | Brain lesion detection | MS lesion tracking, demyelination markers, HLA-DRB1 queries |

### 9.6 Pediatric Imaging Applications

The agent supports pediatric oncology imaging with disease-specific considerations:

- **Posterior Fossa Medulloblastoma:** Pediatric brain tumor imaging in the posterior fossa (most common pediatric brain tumor location). MRI Brain MS Lesion Tracking workflow adapted for pediatric posterior fossa mass detection with age-appropriate measurement thresholds. Triggers genomic queries for medulloblastoma molecular subgroups (WNT, SHH, Group 3, Group 4).
- **ALL Staging CT/MRI:** CT and MRI for staging of acute lymphoblastic leukemia -- hepatosplenomegaly assessment, lymphadenopathy detection, CNS involvement evaluation with brain/spine MRI. Age-appropriate normal organ size references applied.
- **Response Assessment:** Imaging-based treatment response monitoring for pediatric solid tumors using RECIST 1.1 criteria and pediatric-specific adaptations. Longitudinal tracking of tumor volume changes with age-appropriate growth rate thresholds.
- **Sedation Considerations:** Pediatric imaging protocols account for sedation requirements in young children (typically < 6 years for MRI, < 3 years for CT). Protocol recommendations include motion-reduction strategies (rapid sequences, feed-and-wrap for neonates) to minimize sedation needs.

---

## 10. Infrastructure

### 10.1 Technology Stack

| Component | Technology |
|---|---|
| Vector database | Milvus 2.4 (IVF_FLAT, COSINE, 384-dim) |
| Embedding model | BAAI/bge-small-en-v1.5 (sentence-transformers) |
| LLM (primary) | Claude Sonnet 4.6 (Anthropic API) |
| LLM (fallback) | Meta-Llama3-8B-Instruct NIM (OpenAI-compatible) |
| NIM services | VISTA-3D, MAISI, VILA-M3 (NVIDIA Cloud NIMs) |
| DICOM server | Orthanc 24.12.1 (DICOMweb + C-STORE, webhook routing) |
| API framework | FastAPI (async, Prometheus metrics) |
| UI framework | Streamlit (9-tab workbench, NVIDIA themed) |
| Orchestration | Docker Compose (full: 13 services, lite: 6) |
| Workflow models | MONAI bundles + torchxrayvision (pretrained weights) |
| Data models | Pydantic BaseSettings + 20+ models |
| Monitoring | Prometheus metrics endpoint (`/metrics`) |
| Hardware target | DGX Spark (GB10 GPU, 128GB unified, $4,699) |

### 10.2 Service Ports

| Port | Service |
|---|---|
| 4242 | Orthanc DICOM (C-STORE SCP) |
| 8042 | Orthanc DICOMweb REST + web viewer |
| 8520 | NIM LLM (Llama-3 8B Instruct) |
| 8524 | FastAPI REST Server |
| 8525 | Streamlit Imaging Workbench |
| 8530 | NIM VISTA-3D |
| 8531 | NIM MAISI |
| 8532 | NIM VILA-M3 |
| 9091 | Milvus metrics |
| 19530 | Milvus gRPC |

### 10.3 Docker Compose

**Full mode** (13 services): orthanc, milvus-etcd, milvus-minio, milvus-standalone, imaging-api, imaging-streamlit, imaging-setup, nim-llm, nim-vista3d, nim-maisi, nim-vilam3. Requires GPU and NGC_API_KEY.

**Lite mode** (6 services): milvus-etcd, milvus-minio, milvus-standalone, imaging-api, imaging-streamlit, imaging-setup. No GPU required; all NIM features use mock responses.

---

## 11. Deployment

### 11.1 Lite Mode (no GPU)

```bash
cp .env.example .env    # set ANTHROPIC_API_KEY
docker compose -f docker-compose.lite.yml up -d
```

### 11.2 Full Mode (with NVIDIA NIMs)

```bash
cp .env.example .env    # set ANTHROPIC_API_KEY and NGC_API_KEY
docker compose up -d
```

### 11.3 ARM64 Compatibility

All containers use multi-arch images. For NIM containers on DGX Spark (Grace CPU, ARM64), append `-dgx-spark` to standard image tags.

---

## 12. File Structure

```
agent/
+-- api/main.py                         # FastAPI server (port 8524)
|   +-- routes/{meta_agent,nim,workflows,reports,events}.py
+-- app/imaging_ui.py                   # Streamlit 9-tab workbench (port 8525)
+-- config/settings.py                  # Pydantic BaseSettings
+-- src/
|   +-- agent.py                        # Agent orchestrator
|   +-- collections.py                  # 10 Milvus collection schemas
|   +-- cross_modal.py                  # Genomics trigger
|   +-- export.py                       # Markdown, JSON, PDF, FHIR R4
|   +-- knowledge.py                    # Domain knowledge graph
|   +-- metrics.py                      # Prometheus metrics
|   +-- models.py                       # 20+ Pydantic models
|   +-- query_expansion.py             # Term expansion maps
|   +-- rag_engine.py                   # Multi-collection RAG engine
|   +-- scheduler.py                    # APScheduler periodic ingest
|   +-- nim/{base,vista3d_client,maisi_client,vilam3_client,llm_client,service_manager}.py
|   +-- workflows/{base,ct_head_hemorrhage,ct_chest_lung_nodule,ct_coronary_angiography,cxr_rapid_findings,mri_brain_ms_lesion,mri_prostate_pirads}.py
|   +-- ingest/{base,literature_parser,clinical_trials_parser,finding_parser,...,dicom_watcher}.py
+-- tests/                              # Unit tests
+-- scripts/                            # Setup, seeding, validation
+-- docs/                               # Architecture guide, NIM guide, project bible, demo guide
+-- data/sample_images/fullres/         # CXR pathology showcase images
+-- flare/                              # NVIDIA FLARE configs
+-- docker-compose.yml                  # Full (13 services)
+-- docker-compose.lite.yml             # Lite (6 services)
```

---

## 13. DGX Compute Progression

| Phase | Hardware | Price | Scope |
|---|---|---|---|
| **1 -- Proof Build** | DGX Spark | $4,699 | All 6 workflows, mock/cloud NIM fallback |
| **2 -- Departmental** | 1-2x DGX B200 | $500K-$1M | Full NIM stack, PACS integration |
| **3 -- Multi-Site** | 4-8x DGX B200 | $2M-$4M | NVIDIA FLARE federated learning |
| **4 -- AI Factory** | DGX SuperPOD | $7M-$60M+ | Thousands of concurrent studies |

---

## 14. Phase 2 Roadmap

| Item | Priority |
|---|---|
| DICOM SR output (TID 1500 via highdicom) | High |
| NVIDIA FLARE federated learning | Medium |
| Additional MONAI workflows (cardiac, abdominal) | Medium |
| Full PACS integration (dcm4chee) | Medium |
| Grafana + DCGM monitoring dashboard | Low |
| BI-RADS mammography workflow | Low |
| Multi-GPU inference scaling (DGX B200) | Low |

---

## 15. Credits

- **Adam Jones** -- Architecture and implementation
- Part of the **HCLS AI Factory** platform
- Built on NVIDIA DGX Spark ($4,699)
- **Apache 2.0 License**
