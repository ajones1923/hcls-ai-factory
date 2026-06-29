# Precision Autoimmune Intelligence Agent


![Architecture Infographic](docs/images/infographic.jpg)

*Source: [precision-autoimmune-agent](https://github.com/ajones1923/precision-autoimmune-agent)*

**Part of the HCLS AI Factory -- Patient DNA to Drug Candidates**

A multi-collection RAG-powered clinical decision-support agent for autoimmune disease analysis. The agent interprets autoantibody panels, HLA typing, biomarker trends, and genomic data to provide integrated autoimmune assessments including disease activity scoring, flare prediction, and biologic therapy recommendations with pharmacogenomic context. It is designed to surface diagnostic patterns across fragmented clinical records spanning years of multi-specialist visits, addressing the diagnostic odyssey that autoimmune patients commonly face.

---

## Architecture Overview

The agent operates a retrieval-augmented generation (RAG) pipeline built on 14 domain-specific Milvus vector collections. Clinical queries flow through the following stages:

1. **Embedding** -- The input query is encoded using the BGE-small-en-v1.5 model (384 dimensions) with a domain-tuned instruction prefix.

2. **Multi-collection weighted search** -- The query embedding is searched in parallel across all 14 Milvus collections. Each collection carries a configurable relevance weight (e.g., clinical documents at 0.18, patient labs at 0.14, autoantibody panels at 0.12). Results are deduplicated by content hash and capped at 30 evidence items.

3. **Knowledge augmentation** -- The query is matched against a built-in knowledge base of HLA-disease associations (50+ alleles with odds ratios), autoantibody-disease mappings (14+ antibodies with sensitivity/specificity), biologic therapy profiles, and flare biomarker patterns. Matching knowledge items are injected as additional context.

4. **Disease area detection** -- Keyword-based detection identifies relevant disease areas (13 disease categories plus the POTS/hEDS/MCAS triad) to focus retrieval and knowledge augmentation.

5. **Claude synthesis** -- Retrieved evidence, knowledge context, patient context, and conversation history are assembled into a structured prompt and sent to Claude (Sonnet 4) for synthesis. The system prompt enforces citation formats, clinical alert flagging, and cross-domain integration.

6. **Conversation memory** -- A sliding window of recent exchanges (default 3 turns) provides continuity for follow-up questions.

The pipeline supports both synchronous and streaming (SSE) response modes.

### Milvus Collections (14)

| Collection | Weight | Description |
|---|---|---|
| `autoimmune_clinical_documents` | 0.18 | Progress notes, discharge summaries, specialist reports |
| `autoimmune_patient_labs` | 0.14 | Lab results with flags and reference ranges |
| `autoimmune_autoantibody_panels` | 0.12 | ANA, anti-dsDNA, anti-CCP, RF, and 14+ autoantibody types |
| `autoimmune_hla_associations` | 0.08 | HLA allele-disease risk associations with odds ratios |
| `autoimmune_disease_criteria` | 0.08 | ACR/EULAR classification criteria (2010 RA, 2019 SLE, ASAS axSpA) |
| `autoimmune_disease_activity` | 0.07 | DAS28-CRP, SLEDAI-2K, CDAI, BASDAI scoring data |
| `autoimmune_flare_patterns` | 0.06 | Biomarker patterns preceding disease flares |
| `autoimmune_biologic_therapies` | 0.06 | TNF inhibitors, anti-CD20, IL-6R, IL-17A, JAK inhibitors |
| `autoimmune_pgx_rules` | 0.04 | Pharmacogenomic rules (CYP2C19, FCGR3A, HLA-based response) |
| `autoimmune_clinical_trials` | 0.05 | Active and recent clinical trials |
| `autoimmune_literature` | 0.05 | Published literature with year filtering |
| `autoimmune_patient_timelines` | 0.03 | Longitudinal patient event timelines |
| `autoimmune_cross_disease` | 0.02 | Overlap syndromes and cross-disease mechanisms |
| `genomic_evidence` | 0.02 | Shared genomic evidence (read-only, from genomics pipeline) |

---

## Supported Diseases

The agent covers 13 autoimmune diseases:

1. Rheumatoid Arthritis
2. Systemic Lupus Erythematosus (SLE)
3. Multiple Sclerosis
4. Type 1 Diabetes
5. Inflammatory Bowel Disease (Crohn's / Ulcerative Colitis)
6. Psoriasis / Psoriatic Arthritis
7. Ankylosing Spondylitis (Axial Spondyloarthritis)
8. Sjogren's Syndrome
9. Systemic Sclerosis (Scleroderma)
10. Myasthenia Gravis
11. Celiac Disease
12. Graves' Disease
13. Hashimoto's Thyroiditis

The agent also detects the POTS/hEDS/MCAS triad and overlap syndromes (MCTD, rhupus).

---

## Key Features

- **Autoantibody panel interpretation** -- Maps positive autoantibodies to disease associations with sensitivity and specificity data. Supports ANA (with staining patterns), anti-dsDNA, anti-CCP, RF, anti-SSA/SSB, anti-Scl-70, anti-Jo-1, AChR, anti-tTG, and more.

- **HLA association analysis** -- Evaluates HLA typing results against 50+ allele-disease associations with published odds ratios (e.g., HLA-B*27:05 and AS at OR=87.4, HLA-DRB1*04:01 and RA at OR=4.2).

- **Disease activity scoring** -- Calculates standardized scores: DAS28-CRP/ESR for RA, SLEDAI-2K for SLE, CDAI for RA, BASDAI for AS. Reports remission, low, moderate, high, and very high activity levels.

- **Flare prediction** -- Analyzes biomarker patterns (CRP, ESR, IL-6, complement C3/C4, calprotectin, albumin) to predict flare risk on a 0-1 scale with imminent, high, moderate, and low risk levels. Includes contributing and protective factor analysis.

- **Biologic therapy recommendations with PGx context** -- Recommends biologics filtered by diagnosed conditions and pharmacogenomic considerations. Covers TNF inhibitors, anti-CD20, IL-6R, IL-17A, BLyS, IL-12/23, T-cell co-stimulation, and JAK inhibitors with contraindications and monitoring requirements.

- **Diagnostic odyssey analysis** -- Identifies diagnostic patterns across fragmented clinical records spanning years of multi-specialist visits. Detects signals that were individually dismissed but collectively point to a unified autoimmune diagnosis.

- **Overlap syndrome detection** -- Recognizes mixed connective tissue disease (MCTD), rhupus, POTS/hEDS/MCAS triad, and shared pathogenic mechanisms across disease boundaries.

- **FHIR R4 export** -- Exports analysis results in FHIR R4 format for interoperability with electronic health record systems. Also supports Markdown and PDF export.

- **Cross-agent integration** -- Integration points for the Biomarker Agent (inflammation monitoring, longitudinal trends) and Imaging Agent (joint/organ assessment). Event-driven architecture for publishing diagnosis events to other agents.

---

## Quick Start

### Prerequisites

- Python 3.10+
- Running Milvus instance (default: localhost:19530)
- Anthropic API key for Claude

### Docker Compose (recommended)

```bash
# Set your API key
export ANTHROPIC_API_KEY=sk-ant-...

# Start all services (requires shared hcls-network)
docker-compose up -d

# Services:
#   autoimmune-streamlit  -> http://localhost:8531
#   autoimmune-api        -> http://localhost:8532
#   autoimmune-setup      -> one-shot collection creation + seeding
```

### Manual (run.sh)

```bash
# Copy and configure environment
cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY

# Start Streamlit UI only
./run.sh

# Start FastAPI only
./run.sh --api

# Start both UI and API
./run.sh --both

# Create collections and seed knowledge
./run.sh --setup
```

The `run.sh` script automatically creates a Python virtual environment on first run.

### Direct Python

```bash
# Install dependencies
pip install -r requirements.txt

# Start API server
uvicorn api.main:app --host 0.0.0.0 --port 8532 --workers 2

# Start Streamlit UI
streamlit run app/autoimmune_ui.py --server.port 8531 --server.address 0.0.0.0

# Setup collections
python scripts/setup_collections.py --seed
```

---

## Configuration Reference

All settings use the `AUTO_` environment variable prefix and can be set via environment variables or a `.env` file.

| Variable | Default | Description |
|---|---|---|
| `ANTHROPIC_API_KEY` | (none) | Anthropic API key for Claude LLM (required for synthesis) |
| `AUTO_MILVUS_HOST` | `localhost` | Milvus server hostname |
| `AUTO_MILVUS_PORT` | `19530` | Milvus server port |
| `AUTO_STREAMLIT_PORT` | `8531` | Streamlit UI port |
| `AUTO_API_PORT` | `8532` | FastAPI server port |
| `AUTO_LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model identifier |
| `AUTO_LLM_MAX_TOKENS` | `4096` | Maximum tokens in LLM response |
| `AUTO_LLM_TEMPERATURE` | `0.2` | LLM sampling temperature |
| `AUTO_API_KEY` | (empty) | API key for FastAPI authentication; empty disables auth |
| `AUTO_CORS_ORIGINS` | `http://localhost:8080,http://localhost:8531` | Comma-separated CORS origins |
| `AUTO_MAX_REQUEST_SIZE_MB` | `50` | Maximum upload size in MB |
| `AUTO_EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | Sentence transformer model for embeddings |
| `AUTO_EMBEDDING_DIM` | `384` | Embedding vector dimensionality |
| `AUTO_TOP_K_PER_COLLECTION` | `5` | Number of results per collection in search |
| `AUTO_SCORE_THRESHOLD` | `0.40` | Minimum similarity score for search results |
| `AUTO_MAX_EVIDENCE_ITEMS` | `30` | Maximum total evidence items across all collections |
| `AUTO_CONVERSATION_MEMORY_SIZE` | `3` | Number of conversation turns retained |
| `AUTO_MAX_CHUNK_SIZE` | `2500` | Maximum characters per document chunk |
| `AUTO_CHUNK_OVERLAP` | `200` | Character overlap between adjacent chunks |
| `AUTO_CITATION_HIGH` | `0.80` | Score threshold for high-relevance citation |
| `AUTO_CITATION_MEDIUM` | `0.60` | Score threshold for medium-relevance citation |
| `AUTO_FLARE_RISK_IMMINENT` | `0.8` | Flare risk score threshold for imminent level |
| `AUTO_FLARE_RISK_HIGH` | `0.6` | Flare risk score threshold for high level |
| `AUTO_FLARE_RISK_MODERATE` | `0.4` | Flare risk score threshold for moderate level |
| `AUTO_STREAMING_ENABLED` | `true` | Enable streaming responses |
| `AUTO_REQUEST_TIMEOUT_SECONDS` | `60` | Request timeout |
| `AUTO_MILVUS_TIMEOUT_SECONDS` | `10` | Milvus operation timeout |
| `AUTO_LLM_MAX_RETRIES` | `3` | LLM call retry count |
| `AUTO_METRICS_ENABLED` | `true` | Enable Prometheus metrics endpoint |

Collection weights (must sum to approximately 1.0) can also be overridden:

| Variable | Default |
|---|---|
| `AUTO_WEIGHT_CLINICAL_DOCUMENTS` | `0.18` |
| `AUTO_WEIGHT_PATIENT_LABS` | `0.14` |
| `AUTO_WEIGHT_AUTOANTIBODY_PANELS` | `0.12` |
| `AUTO_WEIGHT_HLA_ASSOCIATIONS` | `0.08` |
| `AUTO_WEIGHT_DISEASE_CRITERIA` | `0.08` |
| `AUTO_WEIGHT_DISEASE_ACTIVITY` | `0.07` |
| `AUTO_WEIGHT_FLARE_PATTERNS` | `0.06` |
| `AUTO_WEIGHT_BIOLOGIC_THERAPIES` | `0.06` |
| `AUTO_WEIGHT_PGX_RULES` | `0.04` |
| `AUTO_WEIGHT_CLINICAL_TRIALS` | `0.05` |
| `AUTO_WEIGHT_LITERATURE` | `0.05` |
| `AUTO_WEIGHT_PATIENT_TIMELINES` | `0.03` |
| `AUTO_WEIGHT_CROSS_DISEASE` | `0.02` |
| `AUTO_WEIGHT_GENOMIC_EVIDENCE` | `0.02` |

---

## Demo Patients

The agent ships with 5 synthetic patients, each representing a multi-year autoimmune diagnostic odyssey with realistic clinical documents (progress notes, lab reports, imaging, pathology, genetic testing, referral letters, medication reconciliations, and dismissal documentation).

| Patient | Age/Sex | Primary Condition | Diagnostic Odyssey | Key Pattern |
|---|---|---|---|---|
| Sarah Mitchell | 34F | SLE / Lupus Nephritis | 3.5 years (2022-2025) | Scattered ANA, rising anti-dsDNA, falling C3/C4, proteinuria across 5 specialists |
| Maya Rodriguez | 28F | POTS / hEDS / MCAS Triad | 4 years (2021-2025) | Orthostatic HR spikes in ER vitals, hypermobility never scored, elevated tryptase dismissed |
| David Park | 45M | Ankylosing Spondylitis | 6 years (2019-2025) | Chronic back pain attributed to mechanical causes, HLA-B27 not tested until year 5 |
| Linda Chen | varies | Sjogren's Syndrome | multi-year | Sicca symptoms, anti-SSA/SSB, progressive glandular involvement |
| Rachel Thompson | varies | Mixed Connective Tissue Disease | multi-year | Overlapping features across multiple autoimmune diseases |

Generate demo data with:

```bash
python scripts/generate_demo_patients.py
```

Ingest into Milvus via the API:

```bash
curl -X POST http://localhost:8532/ingest/demo-data
```

---

## Service Ports

| Service | Port | Description |
|---|---|---|
| Streamlit UI | 8531 | Interactive clinical analysis interface |
| FastAPI | 8532 | REST API for programmatic access |

---

## Tech Stack

| Component | Technology |
|---|---|
| Language | Python 3.10 |
| Web framework | FastAPI (API), Streamlit (UI) |
| ASGI server | Uvicorn |
| Vector database | Milvus (pymilvus) |
| Embedding model | BGE-small-en-v1.5 (sentence-transformers) |
| LLM | Claude Sonnet 4 (Anthropic) |
| Data validation | Pydantic v2, pydantic-settings |
| PDF processing | PyPDF2 (reading), ReportLab (generation) |
| Visualization | Plotly, Pandas |
| Monitoring | Prometheus-compatible metrics endpoint |
| Logging | Loguru |
| HTTP client | httpx, requests |
| Container | Docker (multi-stage build, tini init) |
| Orchestration | Docker Compose |

---

## Project Structure

```
precision_autoimmune_agent/
|-- api/
|   |-- main.py                  # FastAPI application and all REST endpoints
|   |-- routes/                  # Route modules (extensible)
|-- app/
|   |-- autoimmune_ui.py         # Streamlit interactive UI
|-- config/
|   |-- settings.py              # Pydantic settings with AUTO_ env prefix
|   |-- logging.py               # Logging configuration
|-- src/
|   |-- agent.py                 # Core AutoimmuneAgent orchestrator
|   |-- models.py                # Pydantic data models (diseases, panels, scores)
|   |-- rag_engine.py            # Multi-collection RAG pipeline
|   |-- collections.py           # Milvus collection manager
|   |-- knowledge.py             # Built-in knowledge base (HLA, antibodies, therapies)
|   |-- diagnostic_engine.py     # Differential diagnosis engine
|   |-- document_processor.py    # PDF ingestion and chunking
|   |-- export.py                # Report export (Markdown, FHIR R4, PDF)
|   |-- timeline_builder.py      # Patient timeline construction
|-- scripts/
|   |-- setup_collections.py     # Collection creation and seeding
|   |-- generate_demo_patients.py# Synthetic patient data generator
|   |-- pdf_engine.py            # PDF document generator for demo data
|   |-- patient_sarah.py         # Sarah Mitchell (SLE) documents
|   |-- patient_maya.py          # Maya Rodriguez (POTS/hEDS/MCAS) documents
|   |-- patients_345.py          # David, Linda, Rachel documents
|   |-- patients_*_enhanced.py   # Enhanced timeline documents
|   |-- patients_dismissals.py   # Diagnostic dismissal documents
|   |-- patients_med_lists.py    # Medication reconciliation documents
|   |-- patients_referrals.py    # Cross-specialist referral letters
|   |-- patients_additional.py   # Additional specialist reports
|-- tests/
|   |-- test_autoimmune.py       # Unit tests
|   |-- test_export.py           # Export tests
|   |-- test_production_readiness.py # Production readiness tests
|-- demo_data/                   # Generated synthetic clinical PDFs
|   |-- sarah_mitchell/
|   |-- maya_rodriguez/
|   |-- david_park/
|   |-- linda_chen/
|   |-- rachel_thompson/
|-- data/                        # Runtime data, cache, reference files
|-- docs/                        # API reference and documentation
|-- docker-compose.yml           # Docker Compose for all services
|-- Dockerfile                   # Multi-stage Docker build
|-- run.sh                       # Startup script (UI, API, both, setup)
|-- requirements.txt             # Python dependencies
|-- .env.example                 # Environment variable template
|-- .streamlit/                  # Streamlit configuration
```

---

## License

This project is licensed under the [Apache License 2.0](LICENSE).

## Author

Adam Jones -- March 2026
