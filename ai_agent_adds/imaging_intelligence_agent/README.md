# Imaging Intelligence Agent

Medical imaging AI agent with RAG knowledge system, NVIDIA NIM integration, and clinical inference workflows. Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory).

## Overview

The Imaging Intelligence Agent provides clinical decision support for radiology through a multi-collection RAG engine backed by 10 Milvus vector collections, 4 NVIDIA NIM microservices, and 4 reference clinical workflows. It searches across medical imaging literature, clinical trials, FDA-cleared devices, acquisition protocols, anatomical references, benchmarks, clinical guidelines, report templates, and public datasets -- synthesizing cross-domain answers grounded in evidence. The agent supports FHIR R4 interoperability for standards-based clinical report exchange, Orthanc DICOM auto-ingestion for automated study processing, cross-modal genomics enrichment linking imaging findings to relevant genetic variants, and NVIDIA FLARE federated learning for privacy-preserving multi-site model training.

### Architecture

```
                    +------------------------------------------------------+
                    |             Streamlit Chat UI (8525)                  |
                    |   Query Input | Evidence Panel | Workflow Runner      |
                    +----------------------------+-------------------------+
                                                 |
                    +----------------------------v-------------------------+
                    |            FastAPI REST Server (8524)                 |
                    |   /api/ask | /query | /search | /reports/generate    |
                    |   /nim/* | /workflows | /workflow/{name}/run         |
                    +--------+------------------+----------------+--------+
                             |                  |                |
              +--------------v----+   +---------v--------+  +---v--------------+
              | RAG Engine        |   | NIM Service Mgr  |  | Workflow Engine   |
              | Multi-collection  |   | Health + Mock    |  | Preprocess->Infer |
              | search + Claude   |   | fallback         |  | ->Postprocess     |
              +--------+----------+   +----+----+---+----+  +--+--+--+--+------+
                       |                   |    |   |          |  |  |  |
              +--------v----------+   +----v-+  | +-v----+   +-v--v--v--v------+
              | 10 Milvus         |   |VISTA |  | |MAISI |   | CT Head | CXR   |
              | Collections       |   |3D    |  | |      |   | CT Lung | MRI   |
              | + genomic_evidence|   |8530  |  | |8531  |   | Hemorrh | Brain |
              | (384-dim BGE)     |   +------+  | +------+   +---------+-------+
              +-------------------+        +----v----+
                                           | VILA-M3 |
              +-------------------+        | 8532    |
              | Knowledge Graph   |        +---------+
              | Pathologies (15)  |
              | Modalities (8)    |   +----v----+
              | Anatomy (15)      |   | LLM NIM |
              +-------------------+   | 8520    |
                                      +---------+
```

## Prerequisites

- Python 3.10+
- Docker and Docker Compose
- NVIDIA GPU (optional -- required only for NIM containers in Full Mode)
- NGC API key (optional -- required only for pulling NIM containers)
- Anthropic API key (for Claude LLM synthesis)

## Quick Start

### Lite Mode (RAG-only, no GPU required)

Runs the RAG knowledge system with Milvus, Streamlit UI, and FastAPI -- all NIM-dependent features use mock responses.

```bash
cp .env.example .env
# Edit .env: set ANTHROPIC_API_KEY

docker compose -f docker-compose.lite.yml up -d
docker compose -f docker-compose.lite.yml logs -f imaging-setup
```

### Full Mode (with NVIDIA NIMs)

Runs the complete stack including VISTA-3D, MAISI, VILA-M3, and Llama-3 NIM containers on GPU.

```bash
cp .env.example .env
# Edit .env: set ANTHROPIC_API_KEY and NGC_API_KEY

docker compose up -d
docker compose logs -f imaging-setup
```

### Service Ports

| Service | Port | Description |
|---|---|---|
| FastAPI REST Server | 8524 | RAG queries, NIM proxy, workflow execution, reports |
| Streamlit Chat UI | 8525 | Interactive chat interface with evidence panel |
| NIM LLM (Llama-3 8B) | 8520 | Text generation and clinical reasoning |
| NIM VISTA-3D | 8530 | 3D medical image segmentation (132 classes) |
| NIM MAISI | 8531 | Synthetic CT volume generation |
| NIM VILA-M3 | 8532 | Vision-language model for radiology |
| Milvus | 19530 | Vector database (gRPC) |
| Orthanc DICOM Server | 8042 / 4242 | DICOM HTTP API (8042) and DICOM C-STORE receiver (4242) |
| Milvus Metrics | 9091 | Health and metrics endpoint |

## Local Development

```bash
# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Ensure Milvus is running (via Docker or standalone)
# Default: localhost:19530

# Create collections and seed reference data
python scripts/setup_collections.py --drop-existing
python scripts/seed_findings.py
python scripts/seed_protocols.py
python scripts/seed_devices.py
python scripts/seed_anatomy.py
python scripts/seed_benchmarks.py
python scripts/seed_guidelines.py
python scripts/seed_report_templates.py
python scripts/seed_datasets.py

# Run the FastAPI server
uvicorn api.main:app --host 0.0.0.0 --port 8524 --reload

# Run the Streamlit UI (in a separate terminal)
streamlit run app/imaging_ui.py --server.port 8525
```

## API Endpoints

### Core RAG

| Method | Path | Description |
|---|---|---|
| POST | `/api/ask` | Meta-agent: question answering with sources and follow-ups |
| POST | `/query` | Full RAG: multi-collection retrieval + LLM synthesis |
| POST | `/search` | Evidence-only search (no LLM) |
| POST | `/find-related` | Cross-collection entity linking |
| GET | `/collections` | List all collections with record counts |
| GET | `/knowledge/stats` | Domain knowledge graph statistics |

### NIM Services

| Method | Path | Description |
|---|---|---|
| GET | `/nim/status` | Health check all NIM services |
| POST | `/nim/vista3d/segment` | Run VISTA-3D 3D segmentation |
| POST | `/nim/maisi/generate` | Generate synthetic CT volume |
| POST | `/nim/vilam3/analyze` | Analyze medical image with VLM |

### Workflows

| Method | Path | Description |
|---|---|---|
| GET | `/workflows` | List all available imaging workflows |
| GET | `/workflow/{name}/info` | Get metadata for a specific workflow |
| POST | `/workflow/{name}/run` | Execute a workflow (mock or live) |

### DICOM Events

| Method | Path | Description |
|---|---|---|
| POST | `/events/dicom-webhook` | Process DICOM study events (auto-triggers workflow) |
| GET | `/events/history` | Paginated DICOM event history |
| GET | `/events/status` | Event bus configuration and routing table |

### Reports

| Method | Path | Description |
|---|---|---|
| POST | `/reports/generate` | Generate clinical report (markdown, JSON, or PDF) |

Export supports Markdown, JSON, PDF, and FHIR R4 DiagnosticReport Bundle.

### Infrastructure

| Method | Path | Description |
|---|---|---|
| GET | `/health` | Service health with collection stats and NIM status |
| GET | `/metrics` | Prometheus-compatible metrics |

## NIM Services

The agent integrates four NVIDIA NIM microservices for on-device medical imaging inference.

| NIM Service | Purpose | Key Capability |
|---|---|---|
| **VISTA-3D** | 3D medical image segmentation | 132 anatomical classes, zero-shot and interactive prompting |
| **MAISI** | Medical AI synthetic imaging | High-resolution synthetic CT with paired segmentation masks |
| **VILA-M3** | Vision-language model | Natural language understanding of radiology images |
| **Llama-3 8B** | Text generation | Clinical reasoning, RAG synthesis, report generation |

All NIM clients inherit from `BaseNIMClient` which provides cached health checks, exponential-backoff retry via tenacity, and automatic mock fallback when services are unavailable. See `docs/NIM_INTEGRATION_GUIDE.md` for setup details.

## Reference Workflows

| Workflow | Modality | Body Region | Target Latency | Model |
|---|---|---|---|---|
| CT Head Hemorrhage Triage | CT | Head | < 90 sec | SegResNet (MONAI wholeBody_ct_segmentation) |
| CT Chest Lung Nodule Tracking | CT | Chest | < 5 min | RetinaNet + SegResNet (MONAI lung_nodule_ct_detection) |
| CXR Rapid Findings | X-ray | Chest | < 30 sec | DenseNet-121 (torchxrayvision CheXpert) |
| MRI Brain MS Lesion Tracking | MRI | Brain | < 5 min | UNEST (MONAI wholeBrainSeg_Large_UNEST_segmentation) |

Each workflow follows a `preprocess -> infer -> postprocess` pipeline with full mock mode support for demonstration and testing.

## Knowledge System

10 Milvus collections with BGE-small-en-v1.5 embeddings (384-dim, IVF_FLAT, COSINE):

| Collection | Content |
|---|---|
| `imaging_literature` | PubMed research papers and reviews |
| `imaging_trials` | ClinicalTrials.gov AI-in-imaging records |
| `imaging_findings` | Imaging finding templates and patterns |
| `imaging_protocols` | Acquisition protocols and parameters |
| `imaging_devices` | FDA-cleared AI/ML medical devices |
| `imaging_anatomy` | Anatomical structure references |
| `imaging_benchmarks` | Model performance benchmarks |
| `imaging_guidelines` | Clinical practice guidelines (ACR, RSNA, NCCN) |
| `imaging_report_templates` | Structured radiology report templates |
| `imaging_datasets` | Public imaging datasets (TCIA, PhysioNet) |
| `genomic_evidence` | *(read-only)* Shared from Stage 2 RAG pipeline |

**Current data:** 2,678 PubMed papers, 12 clinical trials, 124 seed reference records, 3,561,170 genomic evidence vectors (3.56M total vectors across 11 collections).

## Cross-Modal Genomics Integration

When clinical workflows detect high-risk findings (Lung-RADS 4A+ lung nodules, urgent CXR findings), the agent automatically queries the genomic_evidence collection (3.5M vectors from Stage 2 RAG pipeline) for relevant genetic variants (EGFR, ALK, ROS1, KRAS). Results are included in the response as cross-modal enrichment.

## Export Formats

| Format | Function | Output |
|---|---|---|
| Markdown | `export_markdown()` | Structured clinical report |
| JSON | `export_json()` | Full Pydantic model serialization |
| PDF | `export_pdf()` | Professional formatted report (ReportLab) |
| FHIR R4 | `export_fhir()` | DiagnosticReport Bundle with SNOMED CT, LOINC, DICOM coding |

## Testing

```bash
python -m pytest tests/ -v  # 539 tests
python scripts/validate_e2e.py --quick  # 9/9 E2E checks
```

## Project Structure

```
agent/
+-- api/
|   +-- main.py                         # FastAPI server (port 8524)
|   +-- routes/
|       +-- meta_agent.py               # /api/ask endpoint
|       +-- nim.py                      # /nim/* NIM proxy endpoints
|       +-- workflows.py               # /workflow/* endpoints
|       +-- reports.py                  # /reports/generate endpoint
|       +-- events.py                  # DICOM event webhook and history
+-- app/
|   +-- imaging_ui.py                   # Streamlit chat UI (port 8525)
+-- config/
|   +-- settings.py                     # Pydantic BaseSettings configuration
+-- src/
|   +-- agent.py                        # Imaging Intelligence Agent orchestrator
|   +-- collections.py                  # 10 Milvus collection schemas + manager
|   +-- cross_modal.py                  # Cross-modal genomics trigger
|   +-- export.py                       # Report export (PDF, markdown)
|   +-- knowledge.py                    # Domain knowledge graph
|   +-- metrics.py                      # Prometheus metrics
|   +-- models.py                       # Pydantic data models (20+ models)
|   +-- query_expansion.py             # Query expansion maps
|   +-- rag_engine.py                   # Multi-collection RAG engine
|   +-- scheduler.py                    # APScheduler for periodic ingest
|   +-- nim/
|   |   +-- base.py                     # BaseNIMClient ABC
|   |   +-- vista3d_client.py           # VISTA-3D segmentation client
|   |   +-- maisi_client.py            # MAISI synthetic CT client
|   |   +-- vilam3_client.py           # VILA-M3 vision-language client
|   |   +-- llm_client.py             # Llama-3 LLM client
|   |   +-- service_manager.py         # NIM service health manager
|   +-- workflows/
|   |   +-- base.py                     # BaseImagingWorkflow ABC
|   |   +-- ct_head_hemorrhage.py      # CT head hemorrhage triage
|   |   +-- ct_chest_lung_nodule.py    # CT chest lung nodule tracking
|   |   +-- cxr_rapid_findings.py      # CXR rapid findings triage
|   |   +-- mri_brain_ms_lesion.py     # MRI brain MS lesion tracking
|   +-- ingest/
|       +-- base.py                     # Base ingest pipeline
|       +-- literature_parser.py       # PubMed NCBI E-utilities ingest
|       +-- clinical_trials_parser.py  # ClinicalTrials.gov API v2 ingest
|       +-- finding_parser.py          # Imaging finding parser
|       +-- protocol_parser.py         # Protocol parser
|       +-- device_parser.py           # FDA device parser
|       +-- anatomy_parser.py          # Anatomy parser
|       +-- benchmark_parser.py        # Benchmark parser
|       +-- guideline_parser.py        # Guideline parser
|       +-- report_template_parser.py  # Report template parser
|       +-- dicom_watcher.py          # Orthanc DICOM polling service
+-- tests/
+-- docs/
+-- data/
|   +-- sample_images/fullres/         # Full-resolution synthetic CXR images
+-- flare/                              # NVIDIA FLARE federated learning configs
+-- docker-compose.yml                  # Full stack (with NIMs)
+-- docker-compose.lite.yml             # Lite stack (RAG-only, no GPU)
+-- requirements.txt
+-- .env.example
```

## License

Apache 2.0

## Credits

- **Adam Jones** -- Architecture and implementation
- Part of the **HCLS AI Factory** platform
- Built on NVIDIA DGX Spark ($4,699)
