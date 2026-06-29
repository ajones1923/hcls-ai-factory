# Oncology Intelligence Agent


![Architecture Infographic](docs/images/infographic.jpg)

*Source: [precision-oncology-agent](https://github.com/ajones1923/precision-oncology-agent)*

**Closed-loop precision oncology CDSS — from VCF to MTB packet**

Part of the [HCLS AI Factory](https://github.com/your-org/hcls-ai-factory) platform.

---

## Overview

The Oncology Intelligence Agent is a RAG-powered clinical decision support system
designed to accelerate molecular tumor board (MTB) workflows. Given a patient's
somatic and germline variant calls (VCF), the agent identifies actionable
variants, matches them against curated oncology knowledge bases, ranks
candidate therapies, surfaces relevant clinical trials, and assembles a
structured MTB-ready report — all within minutes on commodity hardware.

The agent ingests evidence from CIViC, OncoKB, ClinicalTrials.gov, PubMed,
and NCCN/ASCO/ESMO guidelines into 10 purpose-built Milvus vector collections
plus one shared read-only collection from the upstream genomics pipeline. Every
recommendation is grounded in retrieved evidence with tiered citations so
clinicians can trace the reasoning chain.

Built on the same architecture as the HCLS AI Factory's RAG/Chat pipeline, the
agent uses BGE-small-en-v1.5 embeddings (384-dim), IVF_FLAT/COSINE indexing in
Milvus 2.4, and Claude as the reasoning backbone. A Streamlit UI provides
interactive case exploration, while a FastAPI server exposes programmatic
endpoints for integration with EHR and LIMS systems.

The system is designed to run on a single NVIDIA DGX Spark ($4,699) alongside
the rest of the HCLS AI Factory pipeline, maintaining the project's goal of
end-to-end precision medicine — from Patient DNA to Drug Candidates in under
5 hours.

## Features

- **Variant Annotation** — Automatic lookup of somatic/germline variants against CIViC, OncoKB, and ClinVar with evidence-level tiering (Level 1-4)
- **Therapy Ranking** — Multi-factor ranking of approved and investigational therapies based on variant profile, biomarkers, resistance mechanisms, and guideline concordance
- **Clinical Trial Matching** — Real-time matching against ClinicalTrials.gov based on molecular profile, cancer type, and biomarker eligibility criteria
- **Guideline Retrieval** — Contextual retrieval of NCCN, ASCO, and ESMO recommendations relevant to the patient's molecular and clinical profile
- **Resistance Awareness** — Proactive identification of known resistance mechanisms and bypass strategies for proposed therapies
- **Pathway Analysis** — Cross-referencing of identified variants against signaling pathway maps to uncover druggable nodes and cross-talk
- **MTB Report Generation** — One-click generation of structured molecular tumor board packets with tiered evidence citations
- **Cross-Modal Integration** — Seamless handoff to/from the genomics pipeline (Stage 1) and drug discovery pipeline (Stage 3)
- **FHIR Export** — Export case data and recommendations in FHIR R4 format for EHR integration
- **Scheduled Refresh** — Background scheduler for periodic re-ingestion of CIViC, clinical trials, and literature updates
- **Observability** — Prometheus metrics, OpenTelemetry tracing, and structured logging via Loguru

## Architecture

```
                          +---------------------+
                          |   Streamlit UI      |
                          |   (port 8526)       |
                          +----------+----------+
                                     |
                          +----------v----------+
                          |   FastAPI Server    |
                          |   (port 8527)       |
                          +----------+----------+
                                     |
                 +-------------------+-------------------+
                 |                   |                   |
        +--------v-------+  +-------v--------+  +-------v--------+
        | OncoRAGEngine  |  | TherapyRanker  |  | TrialMatcher   |
        | (multi-coll    |  | (weighted      |  | (biomarker     |
        |  search)       |  |  scoring)      |  |  matching)     |
        +--------+-------+  +-------+--------+  +-------+--------+
                 |                   |                   |
                 +-------------------+-------------------+
                                     |
                          +----------v----------+
                          |  Milvus 2.4         |
                          |  (11 collections)   |
                          |  BGE-small-en-v1.5  |
                          |  IVF_FLAT / COSINE  |
                          +---------------------+
                          |  etcd  |   MinIO    |
                          +--------+-----------+
```

## Collections

| # | Collection          | Description                                       | Source              |
|---|---------------------|---------------------------------------------------|---------------------|
| 1 | `onco_literature`   | PubMed / PMC / preprint chunks by cancer type     | PubMed E-utils      |
| 2 | `onco_trials`       | Clinical trial summaries with biomarker criteria  | ClinicalTrials.gov  |
| 3 | `onco_variants`     | Actionable somatic / germline variants            | CIViC, OncoKB       |
| 4 | `onco_biomarkers`   | Predictive / prognostic biomarkers and assays     | CIViC, literature   |
| 5 | `onco_therapies`    | Approved & investigational therapies with MOA     | OncoKB, FDA labels  |
| 6 | `onco_pathways`     | Signaling pathways, cross-talk, druggable nodes   | KEGG, Reactome      |
| 7 | `onco_guidelines`   | NCCN / ASCO / ESMO guideline recommendations     | Guideline PDFs      |
| 8 | `onco_resistance`   | Resistance mechanisms and bypass strategies       | CIViC, literature   |
| 9 | `onco_outcomes`     | Real-world treatment outcome records              | De-identified RWD   |
|10 | `onco_cases`        | De-identified patient case snapshots              | Synthetic / RWD     |
|11 | `genomic_evidence`  | Read-only VCF-derived evidence (shared)           | Stage 1 pipeline    |

## Quick Start

### Prerequisites

- Docker and Docker Compose v2+
- NVIDIA GPU with CUDA 12+ (recommended for embedding model)
- Anthropic API key

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/hcls-ai-factory.git
cd hcls-ai-factory/ai_agent_adds/precision_oncology_agent/agent

# Configure environment
cp .env.example .env
# Edit .env and set ANTHROPIC_API_KEY
```

### Setup and Run

```bash
# Start all services (Milvus + UI + API + seed)
docker compose up -d

# Watch the one-shot setup / seeding progress
docker compose logs -f onco-setup

# Verify services are healthy
docker compose ps
```

### Access

| Service         | URL                          |
|-----------------|------------------------------|
| Streamlit UI    | http://localhost:8526         |
| FastAPI Docs    | http://localhost:8527/docs    |
| Milvus gRPC     | localhost:19530              |
| Milvus Health   | http://localhost:9091/healthz |

## API Endpoints

| Method | Path                        | Description                                  |
|--------|-----------------------------|----------------------------------------------|
| GET    | `/health`                   | Service health check                         |
| GET    | `/collections/stats`        | Collection entity counts and field lists     |
| POST   | `/query`                    | Free-text oncology RAG query                 |
| POST   | `/cases`                    | Create a new patient case                    |
| GET    | `/cases/{case_id}`          | Retrieve a patient case by ID                |
| POST   | `/cases/{case_id}/analyze`  | Run full MTB analysis on a case              |
| POST   | `/trials/match`             | Match a molecular profile to clinical trials |
| POST   | `/reports/mtb`              | Generate a molecular tumor board packet      |
| POST   | `/reports/fhir`             | Export case data in FHIR R4 format           |
| GET    | `/events/stream`            | SSE stream for real-time analysis updates    |

## Docker Compose Usage

```bash
# Start all services
docker compose up -d

# Start only Milvus (for local development)
docker compose up -d milvus-etcd milvus-minio milvus-standalone

# Run setup/seeding independently
docker compose run --rm onco-setup

# Start only the API server
docker compose up -d onco-api

# View logs
docker compose logs -f onco-streamlit
docker compose logs -f onco-api

# Stop all services
docker compose down

# Stop and remove volumes (full reset)
docker compose down -v
```

## Directory Structure

```
precision_oncology_agent/
  agent/
    api/                    # FastAPI application and route modules
      main.py               #   Application entry point and lifespan
      routes/               #   Route modules (cases, trials, reports, events)
    app/                    # Streamlit UI
      oncology_ui.py        #   Main Streamlit application
    config/                 # Configuration
      settings.py           #   Pydantic settings with ONCO_ env prefix
    data/                   # Runtime data (cache, reference files)
    scripts/                # Setup and seed scripts
    src/                    # Core library
      agent.py              #   Orchestrator / intelligence agent
      case_manager.py       #   Patient case lifecycle
      collections.py        #   Milvus collection schemas and manager
      cross_modal.py        #   Cross-modal pipeline triggers
      export.py             #   FHIR and PDF export
      knowledge.py          #   Knowledge graph utilities
      metrics.py            #   Prometheus metrics
      models.py             #   Pydantic data models
      query_expansion.py    #   Query rewriting and expansion
      rag_engine.py         #   Multi-collection RAG search engine
      scheduler.py          #   APScheduler background tasks
      therapy_ranker.py     #   Multi-factor therapy ranking
      trial_matcher.py      #   Clinical trial matching engine
      ingest/               #   Data ingestion parsers
        civic_parser.py     #     CIViC variant evidence
        oncokb_parser.py    #     OncoKB annotations
        literature_parser.py#     PubMed / PMC literature
        clinical_trials_parser.py  # ClinicalTrials.gov
        guideline_parser.py #     NCCN / ASCO / ESMO guidelines
        pathway_parser.py   #     KEGG / Reactome pathways
        resistance_parser.py#     Resistance mechanism data
        outcome_parser.py   #     Real-world outcome records
      utils/                #   Shared utilities
        vcf_parser.py       #     VCF file parsing (cyvcf2)
        pubmed_client.py    #     PubMed E-utilities client
    tests/                  # Unit and integration tests
    Dockerfile              # Multi-stage Docker build
    docker-compose.yml      # Full stack: Milvus + UI + API + setup
    requirements.txt        # Python dependencies
    LICENSE                 # Apache License 2.0
    README.md               # This file
    .gitignore              # Git ignore rules
    .streamlit/             # Streamlit configuration
      config.toml           #   Server, browser, and theme settings
```

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for the full text.

## Author

**Adam Jones** — 14+ years genomic research experience

Part of the HCLS AI Factory: end-to-end precision medicine platform delivering
Patient DNA to Drug Candidates in under 5 hours on a single NVIDIA DGX Spark.
