# Clinical Imaging Engine -- Project Bible

**Author:** Adam Jones
**Date:** March 2026
**Version:** 3.0.0
**License:** Apache 2.0

---

## 1. Mission

Democratize medical imaging AI by providing a complete, open-source Clinical Imaging Engine (Engine 4) that runs on a single NVIDIA DGX Spark ($4,699). The engine integrates 20 NVIDIA technologies (Community Edition, all free), 9 NIM clients, 9 clinical workflows, 13 Milvus collections (38,028 vectors including 1,938 real PubMed papers), 1,324 tests, AIQ agentic reasoning, NeMo Guardrails, PyRadiomics-CUDA, Holoscan streaming, MONAI Deploy MAPs, and MONAI Label -- combining RAG-grounded clinical knowledge retrieval with on-device NVIDIA NIM inference to support radiology decision-making across CT, MRI, chest X-ray, mammography, ultrasound, and endoscopy modalities. 3-tier deployment: Community/Enterprise/Research.

## 2. Architecture Components

### 2.1 RAG Engine (`src/rag_engine.py`)

Multi-collection retrieval-augmented generation engine:

- Searches 13 imaging-specific Milvus collections in parallel using ThreadPoolExecutor
- Weighted collection scoring (configurable, sums to 1.0)
- Query expansion via domain-specific expansion maps
- Knowledge graph augmentation (pathologies, modalities, anatomy)
- Comparative analysis mode: auto-detects "X vs Y" queries and runs dual retrieval
- LLM synthesis via Claude Sonnet 4.6 (Anthropic API) with Llama-3 NIM fallback
- Conversation memory: injects up to 3 prior exchanges for context

### 2.2 NIM Client Layer (`src/nim/`)

Nine NIM microservice clients with unified interface:

| Client | Service | Port | Capability |
|---|---|---|---|
| `VISTA3DClient` | VISTA-3D | 8530 | 3D segmentation, 132 anatomical classes |
| `MAISIClient` | MAISI | 8531 | Synthetic CT generation with segmentation masks |
| `VILAM3Client` | VILA-M3 | 8532 | Vision-language medical image understanding |
| `LLMClient` | LLM | 8520 | Text generation, clinical reasoning |
| `NVSegmentCTClient` | NV-Segment-CT | -- | CT-specific segmentation |
| `NemotronNanoClient` | Nemotron Nano | -- | Lightweight clinical reasoning |
| `NVGenerateCTClient` | NV-Generate-CT | -- | Synthetic CT volume generation |
| `NVGenerateMRClient` | NV-Generate-MR | -- | Synthetic MR volume generation |
| `NVReasonCXRClient` | NV-Reason-CXR | -- | Chest X-ray reasoning (stub) |

All 9 clients inherit from `BaseNIMClient` (ABC):
- Cached health checks (30-second interval)
- Exponential-backoff retry (3 attempts via tenacity)
- Automatic mock fallback when services are unavailable
- `NIMServiceManager` coordinates all 9 clients and provides aggregated health status

### 2.3 Workflow Engine (`src/workflows/`)

Nine reference clinical workflows following `BaseImagingWorkflow` ABC:
- Pipeline pattern: `preprocess -> infer -> postprocess`
- Full mock mode with clinically realistic synthetic results
- Registered in `WORKFLOW_REGISTRY` for dynamic dispatch

| Workflow | Key | Modality |
|---|---|---|
| CT Head Hemorrhage Triage | `ct_head_hemorrhage` | CT |
| CT Chest Lung Nodule Tracking | `ct_chest_lung_nodule` | CT |
| CT Coronary Angiography | `ct_coronary_angiography` | CT |
| CXR Rapid Findings | `cxr_rapid_findings` | X-ray |
| MRI Brain MS Lesion Tracking | `mri_brain_ms_lesion` | MRI |
| MRI Prostate PI-RADS | `mri_prostate_pirads` | MRI |
| Breast BI-RADS | `breast_birads` | Mammography |
| Thyroid TI-RADS | `thyroid_tirads` | Ultrasound |
| Liver LI-RADS | `liver_lirads` | CT/MRI |

### 2.4 Ingest Pipelines (`src/ingest/`)

Eight domain-specific parsers following `BaseIngestPipeline`:
- Pipeline pattern: `fetch -> parse -> embed -> store`
- Sources: PubMed (NCBI E-utilities), ClinicalTrials.gov (API v2), curated seed data
- Embedding: BGE-small-en-v1.5 (384-dim) via sentence-transformers
- Scheduled via APScheduler (configurable, default weekly)

### 2.5 Knowledge Graph (`src/knowledge.py`)

Domain knowledge providing structured context for RAG augmentation:

| Component | Count | Description |
|---|---|---|
| Pathologies | ~15 | ICD-10 codes, imaging characteristics, severity criteria, AI models |
| Modalities | ~8 | Physics, protocols, typical indications |
| Anatomy | ~15 | VISTA-3D labels, SNOMED codes, FMA IDs, imaging characteristics |

### 2.6 Cross-Modal Genomics Integration (`src/cross_modal.py`)

CrossModalTrigger enriches clinical findings with genomic context via 8 cross-modal triggers:
- Lung-RADS 4A+ findings trigger EGFR, ALK, ROS1, KRAS variant queries
- CXR urgent findings (consolidation) trigger infection genomics queries
- Brain lesion high activity triggers HLA-DRB1 and demyelination marker queries
- CAD-RADS >= 3 triggers LDLR, PCSK9, APOB cardiovascular variant queries
- PI-RADS >= 4 triggers BRCA2, HOXB13 cancer susceptibility queries
- BI-RADS 4+ triggers BRCA1, BRCA2, ATM breast cancer variant queries
- TI-RADS TR4+ triggers RET, BRAF thyroid cancer variant queries
- LI-RADS LR-4+ triggers HFE, SERPINA1 liver disease variant queries
- Queries the shared `genomic_evidence` collection (3.5M vectors)
- Returns genomic hits with variant details, clinical significance, AlphaMissense scores

### 2.7 FHIR R4 Export (`src/export.py`)

Five export formats: Markdown, JSON, PDF, FHIR R4 DiagnosticReport Bundle, and DICOM SR (Structured Report via highdicom TID 1500):
- FHIR Bundle contains Patient, ImagingStudy, Observation, and DiagnosticReport resources
- 54 SNOMED CT codes (103 mapped findings) for findings (hemorrhage=50960005, nodule=416940007, etc.)
- LOINC coding for category (LP29684-5 Radiology) and code (18748-4 Diagnostic imaging study)
- DICOM modality codes (CT, MR, DX, CR, US, PT, MG, RF)
- Observation Interpretation codes mapped from FindingSeverity
- DICOM SR output via highdicom TID 1500 for structured measurement reports

### 2.8 DICOM Ingestion (`src/ingest/dicom_watcher.py`, `api/routes/events.py`)

Orthanc-integrated DICOM auto-ingestion pipeline:
- DicomWatcher polls Orthanc /changes API for StableStudy events
- Webhook endpoint receives study.complete events
- Workflow routing: CT+head->ct_head_hemorrhage, CT+chest->ct_chest_lung_nodule, CT+cardiac->ct_coronary_angiography, CR+chest->cxr_rapid_findings, MR+brain->mri_brain_ms_lesion, MR+prostate->mri_prostate_pirads, MG+breast->breast_birads, US+thyroid->thyroid_tirads, CT/MR+liver->liver_lirads
- Event history with pagination (max 200)

### 2.9 API Layer (`api/`)

FastAPI server on port 8524 with five route modules:

| Module | Prefix | Endpoints |
|---|---|---|
| `meta_agent.py` | `/api` | `/api/ask` -- full question answering with follow-ups |
| `nim.py` | `/nim` | `/nim/status`, `/nim/vista3d/segment`, `/nim/maisi/generate`, `/nim/vilam3/analyze` |
| `workflows.py` | (root) | `/workflows`, `/workflow/{name}/info`, `/workflow/{name}/run` |
| `reports.py` | (root) | `/reports/generate` (markdown, JSON, PDF) |
| `events.py` | `/events` | `/events/dicom-webhook`, `/events/history`, `/events/status` |

Core endpoints on root: `/health`, `/collections`, `/query`, `/search`, `/find-related`, `/knowledge/stats`, `/metrics`, `/demo-cases`, `/demo-cases/{id}/run`, `/protocol/recommend`, `/dose/reference`, `/dose/summary`

### 2.10 Agentic Reasoning (AIQ Toolkit)

AIQ Plan/Execute/Reflect/Refine reasoning with 6 tools for multi-step clinical analysis. The engine plans a series of tool invocations, executes them, reflects on intermediate results, and refines its approach before synthesizing a final answer.

### 2.11 NeMo Guardrails

NeMo Guardrails enforce three safety layers:
- **PII protection** -- Detects and redacts patient identifiers from all outputs
- **Evidence grounding** -- Ensures all claims are traceable to retrieved evidence
- **Disclaimer enforcement** -- Appends clinical disclaimer to all outputs

### 2.12 Radiomics Engine (PyRadiomics-CUDA)

Extracts ~1,500 radiomics features via GPU-accelerated PyRadiomics-CUDA. Features include shape, first-order, GLCM, GLRLM, GLSZM, NGTDM, and GLDM descriptors. Results stored in `imaging_radiomics` collection.

### 2.13 Radiology Report NLP

Full radiology report parsing pipeline: extracts findings, impressions, measurements, and coded diagnoses from free-text radiology reports. Structured results stored in `imaging_reports` collection.

### 2.14 Protocol Optimization

12 ACR indications with patient-specific safety parameters. Recommendations consider age, weight, renal function, contrast allergy history, and pregnancy status. DRL comparison with cumulative dose alerts.

### 2.15 Streaming (Holoscan)

Holoscan real-time streaming pipeline for sub-second AI inference on live ultrasound and endoscopy video feeds.

### 2.16 MONAI Deploy MAPs

9 MONAI Application Packages (MAPs) packaged for clinical deployment following MONAI Deploy containerization standards.

### 2.17 MONAI Label Integration

Interactive annotation with FLARE bridge, enabling radiologists to interactively segment structures and feed corrections back to the model training loop.

### 2.18 Population Analytics (RAPIDS)

GPU-accelerated RAPIDS population analytics for cohort-level imaging trends, disease prevalence monitoring, and outcomes tracking.

### 2.19 UI Layer (`app/imaging_ui.py`)

Streamlit Imaging Workbench on port 8525 with 9 tabs:

| Tab | Description |
|---|---|
| **Evidence Explorer** | Multi-collection RAG queries with pre-filled examples, comparative auto-detection, evidence badges grouped by collection with relevance scores |
| **Workflow Runner** | 9 pre-loaded demo cases (DEMO-001 through DEMO-009), 6-step pipeline animation, annotated AI images, mock/live execution |
| **Image Gallery** | 5 CXR pathology showcase images, cross-modality gallery, 3D volume slice viewer |
| **Protocol Advisor** | Patient-specific imaging protocol recommendations with pre-filled examples, AI dose optimization notes |
| **Device & AI Ecosystem** | 50 FDA-cleared and research AI devices, filterable by modality and clinical task |
| **Dose Intelligence** | 20 protocols with standard vs AI-optimized dose comparison, Plotly charts, 36% average reduction |
| **Reports & Export** | Generate Markdown, JSON, NVIDIA-themed PDF, FHIR R4 DiagnosticReport exports |
| **Patient 360** | Cross-modal genomic enrichment display, interactive Plotly network graph showing gene-finding relationships |
| **Benchmarks & Validation** | Plotly donut chart for collection statistics, model performance benchmarks |

Additional features:
- Sidebar guided tour with 9-step demo flow
- NIM service status indicators
- Multi-turn conversation with context memory
- NVIDIA dark/green themed styling

---

## 3. Milvus Collections

13 imaging-specific collections + 1 read-only cross-agent collection. All use:
- **Index:** IVF_FLAT (nlist=1024, nprobe=16)
- **Metric:** COSINE
- **Embedding dimension:** 384 (BGE-small-en-v1.5)

### Collection Schemas and Seed Counts

| # | Collection | Schema Fields | Seed Source | Seed Count |
|---|---|---|---|---|
| 1 | `imaging_literature` | id, title, text_chunk, source_type, year, modality, body_region, ai_task, disease, keywords, journal | PubMed ingest | 50 |
| 2 | `imaging_trials` | id, title, text_summary, phase, status, sponsor, modality, body_region, ai_task, disease, enrollment, start_year, outcome_summary | ClinicalTrials.gov ingest | 40 |
| 3 | `imaging_findings` | id, text_summary, finding_category, severity, modality, body_region, clinical_significance, differential_diagnosis, recommended_followup, measurement_type, measurement_value, classification_system, classification_score | `seed_findings.py` | 50 |
| 4 | `imaging_protocols` | id, text_summary, protocol_name, modality, body_region, contrast_agent, slice_thickness_mm, radiation_dose, scan_duration, clinical_indication, preprocessing_steps | `seed_protocols.py` | 40 |
| 5 | `imaging_devices` | id, text_summary, device_name, manufacturer, regulatory_status, clearance_date, modality, body_region, ai_task, intended_use, performance_summary, model_architecture | `seed_devices.py` | 50 |
| 6 | `imaging_anatomy` | id, text_summary, structure_name, body_region, system, snomed_code, fma_id, imaging_characteristics, common_pathologies, segmentation_label_id | `seed_anatomy.py` | 30 |
| 7 | `imaging_benchmarks` | id, text_summary, model_name, model_architecture, ai_task, modality, body_region, dataset_name, metric_name, metric_value, training_data_size, inference_time_ms, hardware | `seed_benchmarks.py` | 40 |
| 8 | `imaging_guidelines` | id, text_summary, guideline_name, organization, year, modality, body_region, clinical_indication, classification_system, key_recommendation, evidence_level | `seed_guidelines.py` | 40 |
| 9 | `imaging_report_templates` | id, text_summary, template_name, modality, body_region, finding_type, structured_fields, example_report, coding_system | `seed_report_templates.py` | 50 |
| 10 | `imaging_datasets` | id, text_summary, dataset_name, source, modality, body_region, num_studies, num_images, disease_labels, annotation_type, license_type, download_url | `seed_datasets.py` | 50 |
| 11 | `imaging_radiomics` | id, text_summary, feature_name, feature_class, modality, body_region, extraction_params, clinical_significance, pyradiomics_key | PyRadiomics-CUDA | ~1,500 features |
| 12 | `imaging_reports` | id, text_summary, report_text, findings, impression, measurements, coded_diagnoses, modality, body_region, report_type | Report NLP parser | parsed reports |
| 13 | `imaging_literature` (real) | id, title, text_chunk, source_type, year, modality, body_region, ai_task, disease, keywords, journal, pmid | PubMed ingest (real) | 1,938 real PubMed papers |
| 14 | `genomic_evidence` | *(read-only)* Shared from Stage 2 RAG pipeline | Existing | 3,561,170 |

**Total: 38,028 vectors across 13 imaging collections (including 1,938 real PubMed papers).**

### Collection Search Weights

Configurable in `config/settings.py` (must sum to ~1.0):

| Collection | Weight |
|---|---|
| Literature | 0.18 |
| Findings | 0.15 |
| Trials | 0.12 |
| Guidelines | 0.10 |
| Protocols | 0.08 |
| Devices | 0.08 |
| Benchmarks | 0.08 |
| Anatomy | 0.06 |
| Datasets | 0.05 |
| Report Templates | 0.04 |
| Radiomics | 0.05 |
| Reports | 0.04 |
| Genomic Evidence | 0.03 |

---

## 4. NIM Integration Details

### 4.1 VISTA-3D

- **Image:** `nvcr.io/nvidia/nim/vista3d:latest`
- **Capability:** 3D medical image segmentation supporting 132 anatomical classes
- **Modes:** Zero-shot (automatic), interactive (user-specified classes)
- **Input:** NIfTI volume
- **Output:** Segmentation mask + per-class volumes (mL)
- **GPU memory:** ~8 GB

### 4.2 MAISI

- **Image:** `nvcr.io/nvidia/nim/maisi:latest`
- **Capability:** Medical AI synthetic imaging via latent diffusion
- **Output:** Synthetic CT volume (up to 512x512x512) with paired segmentation masks for up to 127 classes
- **Use cases:** Training data augmentation, algorithm validation
- **GPU memory:** ~12 GB

### 4.3 VILA-M3

- **Image:** `nvcr.io/nvidia/nim/vilam3:latest`
- **Capability:** Vision-language model for medical image understanding
- **Input:** Medical image + natural language question
- **Output:** Natural language findings + confidence score
- **GPU memory:** ~16 GB

### 4.4 Llama-3 8B Instruct

- **Image:** `nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest`
- **Capability:** Text generation, clinical reasoning, RAG synthesis
- **Endpoint:** OpenAI-compatible `/v1/chat/completions`
- **GPU memory:** ~16 GB

### 4.5 Additional NIM Clients (v2.1+)

| NIM Client | Capability |
|---|---|
| NV-Segment-CT | CT-specific segmentation with organ-level granularity |
| Nemotron Nano | Lightweight clinical reasoning for edge deployment |
| NV-Generate-CT | Synthetic CT volume generation |
| NV-Generate-MR | Synthetic MR volume generation |
| NV-Reason-CXR | Chest X-ray reasoning (stub, planned) |

### 4.6 Cloud NIM Fallback

When local NIM containers are unavailable, the engine falls back to NVIDIA Cloud NIMs:
- `meta/llama-3.1-8b-instruct` via `integrate.api.nvidia.com`
- `meta/llama-3.2-11b-vision-instruct` via `integrate.api.nvidia.com`

---

## 5. Reference Workflows

### 5.1 CT Head Hemorrhage Triage

| Attribute | Value |
|---|---|
| Target latency | < 90 seconds |
| Sensitivity target | > 95% for hemorrhage > 5 mL |
| Model | SegResNet (MONAI wholeBody_ct_segmentation) |
| Preprocessing | RAS orientation, 1mm isotropic, blood window (0-80 HU) |
| Postprocessing | Volume estimation, midline shift, max thickness |
| Urgency | BTF thresholds: >30mL OR >5mm shift = Critical; >5mL = Urgent |
| Output | Finding + measurements + WorklistEntry |

### 5.2 CT Chest Lung Nodule Tracking

| Attribute | Value |
|---|---|
| Target latency | < 5 minutes |
| Detection target | > 90% for nodules >= 4 mm |
| Detection model | RetinaNet (MONAI) |
| Segmentation model | SegResNet (MONAI) |
| Longitudinal | SyN diffeomorphic registration, volume doubling time |
| Classification | ACR Lung-RADS v2022 |
| Cross-modal trigger | Lung-RADS 4A+ triggers genomic variant queries (EGFR, ALK, ROS1, KRAS) |

### 5.3 CT Coronary Angiography

| Attribute | Value |
|---|---|
| Target latency | < 3 minutes |
| Measurements | Calcium score (Agatston), stenosis grading per vessel |
| Classification | CAD-RADS |
| Plaque analysis | Low-attenuation plaque, positive remodeling detection |
| Vessel coverage | LAD, LCx, RCA, Left Main |
| Cross-modal trigger | High-risk findings trigger FH gene queries (LDLR, PCSK9, APOB) |

### 5.4 CXR Rapid Findings

| Attribute | Value |
|---|---|
| Target latency | < 30 seconds |
| Model | DenseNet-121 multi-label classification (CheXpert pretrained) |
| Findings | Consolidation, effusion, pneumothorax, cardiomegaly, atelectasis, edema, nodule |
| Localization | GradCAM heatmap overlay |
| Output | Multi-label findings with confidence scores |

### 5.5 MRI Brain MS Lesion Tracking

| Attribute | Value |
|---|---|
| Target latency | < 5 minutes |
| Model | UNEST (MONAI wholeBrainSeg_Large_UNEST_segmentation) |
| Registration | SyN diffeomorphic (ANTsPy) |
| Longitudinal | Lesion matching, new/enlarged/resolved classification |
| Metrics | Total lesion volume, lesion count, lesion change rate |
| Disease activity | Stable / Active / Highly Active |

### 5.6 MRI Prostate PI-RADS

| Attribute | Value |
|---|---|
| Target latency | < 5 minutes |
| Classification | PI-RADS v2.1 scoring |
| Assessment | Lesion characterization, zone localization |
| Output | PI-RADS score, lesion details, biopsy recommendation |

### 5.7 Breast BI-RADS (`breast_birads`)

| Attribute | Value |
|---|---|
| Target latency | < 5 minutes |
| Classification | BI-RADS 0-6 |
| Assessment | Mass characterization, calcification analysis, architectural distortion |
| Cross-modal trigger | BI-RADS 4+ triggers BRCA1, BRCA2, ATM variant queries |

### 5.8 Thyroid TI-RADS (`thyroid_tirads`)

| Attribute | Value |
|---|---|
| Target latency | < 3 minutes |
| Classification | ACR TI-RADS TR1-TR5 |
| Assessment | Nodule composition, echogenicity, shape, margin, echogenic foci |
| Cross-modal trigger | TI-RADS TR4+ triggers RET, BRAF thyroid cancer variant queries |

### 5.9 Liver LI-RADS (`liver_lirads`)

| Attribute | Value |
|---|---|
| Target latency | < 5 minutes |
| Classification | LI-RADS LR-1 to LR-5, LR-M, LR-TIV |
| Assessment | Arterial phase hyperenhancement, washout, capsule, threshold growth |
| Cross-modal trigger | LI-RADS LR-4+ triggers HFE, SERPINA1 liver disease variant queries |

---

## 6. API Endpoints Catalog

### Core RAG Endpoints

| Method | Path | Request Body | Response |
|---|---|---|---|
| POST | `/api/ask` | `AskRequest` (question, modality, body_region, top_k, conversation_history) | `AskResponse` (answer, sources, follow_up_questions) |
| POST | `/query` | `QueryRequest` (question, modality, body_region, top_k, collections, year_min/max) | `QueryResponse` (answer, evidence_count, collections_searched) |
| POST | `/search` | `SearchRequest` (question, modality, top_k, collections) | `SearchResponse` (hits, total_hits, knowledge_context) |
| POST | `/find-related` | `FindRelatedRequest` (entity, top_k) | `FindRelatedResponse` (collections with hits) |

### NIM Endpoints

| Method | Path | Request Body | Response |
|---|---|---|---|
| GET | `/nim/status` | -- | `NIMStatusResponse` (services, available/mock/unavailable counts) |
| POST | `/nim/vista3d/segment` | `SegmentRequest` (input_path, target_classes) | `SegmentResponse` (classes, volumes, inference_time) |
| POST | `/nim/maisi/generate` | `GenerateRequest` (body_region, resolution, num_classes) | `GenerateResponse` (resolution, generation_time) |
| POST | `/nim/vilam3/analyze` | `AnalyzeRequest` (question, input_path) | `AnalyzeResponse` (answer, findings, confidence) |

### Workflow Endpoints

| Method | Path | Request Body | Response |
|---|---|---|---|
| GET | `/workflows` | -- | `WorkflowListResponse` (workflow metadata) |
| GET | `/workflow/{name}/info` | -- | `WorkflowInfo` (modality, body_region, target_latency, models) |
| POST | `/workflow/{name}/run` | `WorkflowRunRequest` (input_path, mock_mode) | `WorkflowRunResponse` (findings, measurements, classification, severity) |

### Demo Case Endpoints

| Method | Path | Request Body | Response |
|---|---|---|---|
| GET | `/demo-cases` | -- | List of 9 demo case metadata |
| POST | `/demo-cases/{id}/run` | -- | Workflow result + genomic context + talking points |

### Protocol and Dose Endpoints

| Method | Path | Request Body | Response |
|---|---|---|---|
| POST | `/protocol/recommend` | indication, patient_age, patient_weight_kg | Protocol recommendation with AI optimization |
| GET | `/dose/reference` | -- | All 20 dose comparison protocols |
| GET | `/dose/summary` | -- | Summary statistics (avg/max/min reduction by modality) |
| GET | `/dose/comparison/{region}` | -- | Dose comparison for specific body region |

### Report Endpoints

| Method | Path | Request Body | Response |
|---|---|---|---|
| POST | `/reports/generate` | `ReportRequest` (question, format: markdown/json/pdf) | `ReportResponse` or PDF binary |

### DICOM Event Endpoints

| Method | Path | Request Body | Response |
|---|---|---|---|
| POST | `/events/dicom-webhook` | `DicomStudyEvent` (event_type, study_uid, patient_id, modality, body_region) | `DicomIngestionResult` (workflow_triggered, workflow_result) |
| GET | `/events/history` | Query: limit, offset | Paginated event list |
| GET | `/events/status` | -- | Routing table, history count, Orthanc URL |

### Infrastructure Endpoints

| Method | Path | Response |
|---|---|---|
| GET | `/health` | `HealthResponse` (status, collections, total_vectors, nim_services) |
| GET | `/collections` | `List[CollectionInfo]` (name, count, label) |
| GET | `/knowledge/stats` | Knowledge graph statistics |
| GET | `/metrics` | Prometheus-compatible metrics |

---

## 7. UI Features

### Streamlit Imaging Workbench (port 8525)

**9-tab interface:**

1. **Evidence Explorer** -- Multi-collection RAG queries with pre-filled example queries. Auto-detects comparative queries ("CT vs MRI") and produces structured comparison tables. Evidence sources display collection badges and cosine similarity scores.

2. **Workflow Runner** -- Execute 9 reference workflows via 9 pre-loaded demo cases (DEMO-001 through DEMO-009). Features 6-step pipeline animation showing the processing stages, annotated AI images with overlay results, and mock/live mode toggle.

3. **Image Gallery** -- Showcase of 5 CXR pathology images (cardiomegaly, pleural effusion, pneumothorax, consolidation, atelectasis) with AI annotations. Cross-modality gallery and 3D volume slice viewer for CT/MRI data.

4. **Protocol Advisor** -- Patient-specific imaging protocol recommendations based on clinical indication, age, and weight. Pre-filled example queries. AI dose optimization notes with DLIR recommendations.

5. **Device & AI Ecosystem** -- Searchable catalog of 50 FDA-cleared and research AI devices. Filterable by modality (CT, MRI, X-ray, etc.) and clinical task (detection, segmentation, triage).

6. **Dose Intelligence** -- Dashboard comparing standard vs AI-optimized radiation doses across 20 protocols. Plotly bar charts, 36% average dose reduction, DLIR technique details, ALARA compliance.

7. **Reports & Export** -- Generate clinical reports in 5 formats: Markdown (copy-paste), JSON (programmatic), NVIDIA-themed PDF (documentation), FHIR R4 DiagnosticReport Bundle (EHR interoperability with 54 SNOMED CT codes (103 mapped findings) and LOINC coding), and DICOM SR (Structured Report via highdicom TID 1500).

8. **Patient 360** -- Cross-modal genomic enrichment results showing gene-finding relationships. Interactive Plotly network graph connecting imaging findings to genetic risk factors. Gene details with clinical significance and AlphaMissense scores.

9. **Benchmarks & Validation** -- Plotly donut chart visualizing collection statistics. Model performance benchmarks across workflows. Validation metrics and quality indicators.

**Additional features:**
- Sidebar guided tour with 9-step demo flow
- NIM service status indicators (available / mock / unavailable)
- Multi-turn conversation with context memory (up to 3 exchanges)
- NVIDIA dark/green themed styling throughout

---

## 8. Docker Deployment

### Full Stack (`docker-compose.yml`)

13 services total:
- `milvus-etcd` -- etcd key-value store for Milvus metadata
- `milvus-minio` -- MinIO object storage for Milvus data
- `milvus-standalone` -- Milvus 2.4 vector database
- `imaging-streamlit` -- Streamlit Imaging Workbench (9 tabs)
- `imaging-api` -- FastAPI REST server
- `imaging-setup` -- One-shot collection creation and data seeding (38,028 vectors)
- `orthanc` -- Orthanc DICOM server (PACS, port 8042 HTTP / 4242 DICOM)
- `nim-llm` -- Meta Llama-3 8B Instruct
- `nim-vista3d` -- NVIDIA VISTA-3D segmentation
- `nim-maisi` -- NVIDIA MAISI synthetic imaging
- `nim-vilam3` -- VILA-M3 vision-language model

### Lite Stack (`docker-compose.lite.yml`)

6 services (no GPU required):
- Milvus infrastructure (etcd, MinIO, standalone)
- imaging-streamlit (NIM_MODE=mock)
- imaging-api (NIM_MODE=mock)
- imaging-setup

### Volumes

| Volume | Purpose |
|---|---|
| `etcd_data` | Milvus metadata persistence |
| `minio_data` | Milvus object storage persistence |
| `milvus_data` | Milvus vector data persistence |
| `nim_models` | Cached NIM model weights (full stack only) |
| `orthanc_data` | Orthanc DICOM storage persistence |

---

## 9. Configuration

All configuration via Pydantic `BaseSettings` in `config/settings.py`. Environment variables use the `IMAGING_` prefix. Key settings:

| Setting | Default | Description |
|---|---|---|
| `MILVUS_HOST` | localhost | Milvus server hostname |
| `MILVUS_PORT` | 19530 | Milvus gRPC port |
| `EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Embedding model name |
| `EMBEDDING_DIMENSION` | 384 | Embedding vector dimension |
| `LLM_PROVIDER` | anthropic | LLM provider (anthropic or nim) |
| `ANTHROPIC_API_KEY` | -- | Claude API key (required for primary LLM) |
| `NIM_MODE` | local | NIM mode: local (Docker) or mock |
| `NIM_ALLOW_MOCK_FALLBACK` | true | Fall back to mock when NIM unavailable |
| `TOP_K_PER_COLLECTION` | 5 | Default results per collection |
| `SCORE_THRESHOLD` | 0.4 | Minimum cosine similarity score |
| `INGEST_SCHEDULE_HOURS` | 168 | Periodic ingest interval (weekly) |
| `MAX_CONVERSATION_CONTEXT` | 3 | Prior exchanges for conversation memory |

---

## 10. Data Statistics

| Metric | Count |
|---|---|
| Vectors (13 collections) | 38,028 |
| Real PubMed papers | 1,938 |
| Genomic evidence vectors | 3,561,170 |
| Milvus collections (owned) | 13 |
| Milvus collections (shared) | 1 (genomic_evidence, read-only) |
| Reference workflows | 9 |
| NIM clients | 9 |
| Demo cases | 9 (DEMO-001 through DEMO-009) |
| Tests | 1,324 |
| NVIDIA technologies integrated | 20 (Community Edition, all free) |
| SNOMED codes | 103 |
| Cross-modal triggers | 8 |
| Scoring systems | 7 (Lung-RADS, BI-RADS, TI-RADS, LI-RADS, CAD-RADS, PI-RADS, ASPECTS) |
| Export formats | 5 (Markdown, JSON, PDF, FHIR R4, DICOM SR) |
| MONAI Deploy MAPs | 9 |
| Radiomics features | ~1,500 (PyRadiomics-CUDA) |
| ACR protocol indications | 12 |
| UI tabs | 9 |
| React portal pages | 10 |
| Docker services (full) | 13 |
| Docker services (lite) | 6 |
| CXR pathology images | 5 |
| Dose comparison protocols | 20 |
| AI device catalog entries | 50 |
| Intelligence agents in HCLS AI Factory | 11 |

---

## 11. Phase 2 Roadmap

### Implemented (Phase 1.1 -- March 2026)

| Feature | Description | Status |
|---|---|---|
| FHIR R4 Output | DiagnosticReport Bundle with SNOMED CT, LOINC, DICOM coding | Implemented |
| DICOM Server Integration | Orthanc DICOM server with webhook auto-routing | Implemented |
| Cross-Modal Triggers | Lung-RADS 4A+ triggers genomic queries (EGFR/ALK/ROS1/KRAS) | Implemented |
| Federated Learning | NVIDIA FLARE configs: 3 jobs, multi-site, mTLS, HE (see `flare/`) | Implemented |
| Cloud NIM Inference | NVIDIA Cloud NIMs via integrate.api.nvidia.com (Llama-3.1-8B + Llama-3.2-11B-Vision) | Implemented |
| Real Pretrained Weights | CheXpert DenseNet-121, MONAI RetinaNet, SegResNet, UNEST | Implemented |
| 9-Tab Streamlit Workbench | Evidence Explorer, Workflow Runner, Image Gallery, Protocol Advisor, Device & AI Ecosystem, Dose Intelligence, Reports & Export, Patient 360, Benchmarks & Validation | Implemented |
| Image Gallery | 5 CXR pathology showcase, 3D volume viewer | Implemented |
| Demo Cases | 4 pre-loaded clinical scenarios (DEMO-001 through DEMO-004) | Implemented |
| 6 Workflows | Added CT Coronary Angiography and MRI Prostate PI-RADS | Implemented |
| Guided Tour | Sidebar 9-step demo walkthrough | Implemented |

### Implemented (Phase 2.0 -- April 2026)

| Feature | Description | Status |
|---|---|---|
| DICOM SR Output | Structured reports via highdicom TID 1500 | Implemented |
| Population Analytics | GPU-accelerated RAPIDS cohort-level imaging trends and outcomes | Implemented |
| AIQ Agentic Reasoning | Plan/Execute/Reflect/Refine with 6 tools | Implemented |
| NeMo Guardrails | PII protection, evidence grounding, disclaimer enforcement | Implemented |
| 9 NIM Clients | Added NV-Segment-CT, Nemotron Nano, NV-Generate-CT, NV-Generate-MR, NV-Reason-CXR (stub) | Implemented |
| 9 Workflows | Added Breast BI-RADS, Thyroid TI-RADS, Liver LI-RADS | Implemented |
| 9 Demo Cases | DEMO-005 through DEMO-009 for new workflows | Implemented |
| Radiomics | ~1,500 features via PyRadiomics-CUDA, imaging_radiomics collection | Implemented |
| Report NLP | Full radiology report parsing pipeline, imaging_reports collection | Implemented |
| 38,028 Vectors | 1,938 real PubMed papers, 13 collections total | Implemented |
| 1,324 Tests | Expanded from 620 across all new modules | Implemented |
| 54 SNOMED CT Codes (103 mapped findings) | Expanded from 34 for comprehensive finding coverage | Implemented |
| 8 Cross-Modal Triggers | Added BI-RADS, TI-RADS, LI-RADS triggers | Implemented |
| Protocol Optimization | 12 ACR indications with patient-specific safety | Implemented |
| Dose Tracking | DRL comparison, cumulative alerts | Implemented |
| Holoscan Streaming | Real-time ultrasound/endoscopy pipeline | Implemented |
| MONAI Deploy MAPs | 9 MAPs packaged for clinical deployment | Implemented |
| MONAI Label | Interactive annotation with FLARE bridge | Implemented |
| 3D Visualization | Three.js rotating point cloud | Implemented |
| React Portal | Full React portal with 10 pages | Implemented |
| 20 NVIDIA Technologies | Community Edition, all free | Implemented |
| Live DICOM Analysis | Real GPU inference with DenseNet-121 (CXR) + 6 MONAI bundles (1.87 GB), drag-and-drop upload | Implemented |
| 3-Tier Deployment | Community/Enterprise/Research | Implemented |

### Planned (Phase 3)

| Feature | Description | Status |
|---|---|---|
| Live Local NIM Inference | Local GPU deployment of all 9 NIM services | NIM clients ready |
| Multi-GPU Scaling | DGX B200 multi-GPU inference | Planned |
| Domain-Specific Embeddings | RadBERT fine-tuned embeddings for radiology | Planned |
| Multi-Institutional Validation | Prospective evaluation of cross-modal triggers | Planned |

---

*This document serves as the complete implementation specification for the Clinical Imaging Engine (Engine 4). For architecture details see the Design Document. For demo walkthrough see the Demo Guide.*
