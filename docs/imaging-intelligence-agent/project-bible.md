# Imaging Intelligence Agent — Project Bible

> **Purpose:** Complete implementation reference for building the HCLS Imaging Intelligence Agent on NVIDIA DGX Spark. Import this document into a Claude Code session as context for implementation.
>
> **License:** Apache 2.0 | **Author:** Adam Jones | **Date:** February 2026

---

## Table of Contents

1. [Project Overview & Goals](#1-project-overview--goals)
2. [DGX Spark Hardware Reference](#2-dgx-spark-hardware-reference)
3. [Repository Layout](#3-repository-layout)
4. [Docker Compose Services](#4-docker-compose-services)
5. [PostgreSQL + pgvector Schema](#5-postgresql--pgvector-schema)
6. [Pydantic Data Models](#6-pydantic-data-models)
7. [Orthanc DICOM Server Configuration](#7-orthanc-dicom-server-configuration)
8. [MONAI Deploy MAP Pattern](#8-monai-deploy-map-pattern)
9. [CT Head Hemorrhage Workflow](#9-ct-head-hemorrhage-workflow)
10. [CT Chest Lung Nodule Workflow](#10-ct-chest-lung-nodule-workflow)
11. [CXR Rapid Findings Workflow](#11-cxr-rapid-findings-workflow)
12. [MRI Brain MS Lesion Workflow](#12-mri-brain-ms-lesion-workflow)
13. [DICOM SR Output (highdicom)](#13-dicom-sr-output-highdicom)
14. [FHIR DiagnosticReport Output](#14-fhir-diagnosticreport-output)
15. [LangGraph Agent Architecture](#15-langgraph-agent-architecture)
16. [NIM LLM Integration](#16-nim-llm-integration)
17. [Embedding Service](#17-embedding-service)
18. [Nextflow Pipeline Orchestration](#18-nextflow-pipeline-orchestration)
19. [Streamlit Portal](#19-streamlit-portal)
20. [Monitoring Stack](#20-monitoring-stack)
21. [Testing Strategy](#21-testing-strategy)
22. [ARM64 Compatibility Guide](#22-arm64-compatibility-guide)
23. [Configuration Reference](#23-configuration-reference)
24. [Implementation Sequence](#24-implementation-sequence)

---

## 1. Project Overview & Goals

### What This Agent Does

The Imaging Intelligence Agent processes CT, MRI, and X-ray studies using NVIDIA MONAI models on DGX hardware. It provides automated detection, segmentation, classification, longitudinal tracking, and clinical triage — from DICOM ingestion through agentic inference to structured clinical output.

### Four Reference Workflows

| Workflow | Modality | Target Latency | Key Metric |
|---|---|---|---|
| Hemorrhage Triage | CT Head | < 90 seconds | Sensitivity > 95% for hemorrhage > 5 mL |
| Lung Nodule Tracking | CT Chest | < 5 minutes | Detection > 90% for nodules >= 4 mm |
| Rapid Findings | CXR | < 30 seconds | Pneumothorax sensitivity > 95% |
| MS Lesion Tracking | MRI Brain | < 5 minutes | 3D U-Net on FLAIR sequences |

### Pipeline Pattern

Every workflow follows the same pattern:
1. DICOM ingestion (DICOMweb STOW-RS or DIMSE C-STORE)
2. Event trigger (study.complete)
3. Prior study retrieval (longitudinal comparison)
4. GPU inference (MONAI Deploy MAP on DGX Spark)
5. Post-processing (measurements, classifications, embeddings)
6. Persistence (findings → PostgreSQL, embeddings → pgvector)
7. Output encoding (DICOM SR → PACS, FHIR → EHR)
8. Triage routing (worklist prioritization, alerts)
9. Provenance (immutable audit trail)

### HCLS AI Factory Integration

This agent is one node in the broader HCLS AI Factory. Cross-modal triggers include:
- **Imaging → Genomics (Parabricks):** Lung-RADS 4B+ triggers tumor profiling
- **Imaging → Drug Discovery (BioNeMo):** Quantitative imaging endpoints for trials
- **Imaging → Clinical Reasoning (NIM LLM):** RAG-grounded clinical reports
- **Imaging → Biomarker Agent:** Genomic + imaging biomarker fusion

---

## 2. DGX Spark Hardware Reference

### Specifications

| Parameter | Value |
|---|---|
| CPU | NVIDIA Grace (ARM64 / aarch64) |
| GPU | NVIDIA Blackwell GB10, 1 GPU |
| Memory | 128 GB unified LPDDR5x (CPU + GPU shared pool) |
| Storage | Up to 4 TB NVMe |
| Storage Access | GPUDirect Storage (zero-copy GPU access) |
| Price | $3,999 |
| OS | Ubuntu-based (NVIDIA DGX OS) |
| NVAIE Cost | Zero at desktop-class |

### Critical: ARM64 Architecture

**ALL containers must be ARM64-compatible.** The Grace CPU is aarch64, not x86_64. This affects:
- Base Docker images (must use ARM64 variants)
- Python wheel availability (most scientific packages have ARM64 wheels)
- NIM containers (use `-dgx-spark` variant image tags)
- Any compiled C/C++ extensions

### Unified Memory Model

The 128 GB is **shared** between CPU and GPU — there is no separate GPU VRAM. This means:
- No explicit CPU→GPU data transfers needed for many operations
- Memory pressure from CPU workloads reduces GPU-available memory
- Monitor total system memory, not just "GPU memory"

### GPUDirect Storage

NVMe → GPU memory with zero CPU bounce buffers. Configure with:
```bash
# Verify GPUDirect Storage support
nvidia-smi topo -m
# Check NVMe devices
nvme list
```

### DGX Compute Progression

| Phase | Hardware | Price | Scope |
|---|---|---|---|
| 1 — Proof Build | DGX Spark | $3,999 | 1-2 workflows |
| 2 — Departmental | 1-2x DGX B200 | $500K-$1M | All workflows, PACS integration |
| 3 — Multi-Site | 4-8x DGX B200 + InfiniBand | $2M-$4M | FLARE federated learning |
| 4 — AI Factory | DGX SuperPOD | $7M-$60M+ | Thousands of concurrent studies |

---

## 3. Repository Layout

```
hls-imaging-agent/
├── main.nf                          # Nextflow DSL2 entry point
├── nextflow.config                  # Profiles: docker, dgx_spark
├── docker-compose.yml               # All services
├── docker-compose.dev.yml           # Dev overrides (mock NIM, etc.)
├── start-services.sh                # Service startup script
├── demo.sh                          # Demo launcher
├── .env.example                     # Environment variable template
├── requirements.txt                 # Python dependencies (portal, agent)
│
├── modules/                         # Nextflow workflow modules
│   ├── ct_head_hemorrhage.nf
│   ├── ct_chest_lung_nodule.nf
│   ├── cxr_rapid_findings.nf
│   └── mri_brain_ms_lesion.nf
│
├── maps/                            # MONAI Deploy Application Packages
│   ├── ct_head_hemorrhage/
│   │   ├── Dockerfile
│   │   ├── app.py                   # Application class
│   │   ├── operators.py             # Inference + post-processing operators
│   │   └── requirements.txt
│   ├── ct_chest_lung_nodule/
│   │   ├── Dockerfile
│   │   ├── app.py
│   │   ├── operators.py
│   │   └── requirements.txt
│   ├── cxr_rapid_findings/
│   │   ├── Dockerfile
│   │   ├── app.py
│   │   ├── operators.py
│   │   └── requirements.txt
│   └── mri_brain_ms_lesion/
│       ├── Dockerfile
│       ├── app.py
│       ├── operators.py
│       └── requirements.txt
│
├── agent/                           # LangGraph clinical reasoning agent
│   ├── __init__.py
│   ├── graph.py                     # StateGraph definition
│   ├── nodes.py                     # Agent node implementations
│   ├── tools.py                     # MCP tool definitions
│   ├── state.py                     # AgentState TypedDict
│   └── prompts.py                   # Agent persona system prompts
│
├── services/                        # Microservices
│   ├── dicom_listener/
│   │   ├── Dockerfile
│   │   ├── listener.py              # DICOMweb + Orthanc webhook listener
│   │   └── requirements.txt
│   ├── fhir_publisher/
│   │   ├── Dockerfile
│   │   ├── publisher.py             # FHIR DiagnosticReport output
│   │   └── requirements.txt
│   ├── embedding_service/
│   │   ├── Dockerfile
│   │   ├── embedder.py              # Study/series/lesion embeddings
│   │   └── requirements.txt
│   └── portal/
│       ├── Dockerfile
│       ├── app.py                   # Streamlit dashboard
│       └── requirements.txt
│
├── src/                             # Shared library code
│   ├── __init__.py
│   ├── models.py                    # Pydantic data models
│   ├── db.py                        # PostgreSQL + pgvector client
│   ├── dicom_utils.py               # DICOM SR, GSPS, SEG helpers
│   ├── fhir_utils.py                # FHIR resource builders
│   └── config.py                    # Environment-based configuration
│
├── db/
│   ├── init.sql                     # PostgreSQL + pgvector schema
│   └── migrations/                  # Schema migration scripts
│
├── models/
│   ├── README.md                    # Model download instructions
│   └── download_models.sh           # MONAI Model Zoo download script
│
├── config/
│   ├── orthanc.json                 # Orthanc DICOM server config
│   ├── prometheus.yml               # Monitoring config
│   └── grafana/
│       └── dashboards/
│           └── imaging-agent.json   # Grafana dashboard
│
├── tests/
│   ├── conftest.py                  # Shared fixtures
│   ├── test_ct_head.py
│   ├── test_ct_chest.py
│   ├── test_cxr.py
│   ├── test_mri_brain.py
│   ├── test_db.py
│   ├── test_dicom_sr.py
│   └── test_fhir.py
│
├── scripts/
│   ├── download_models.sh           # MONAI Model Zoo downloads
│   └── seed_test_data.sh            # Synthetic DICOM test data
│
└── docs/
    └── diagrams/
        └── architecture.mmd         # Mermaid architecture diagram
```

---

## 4. Docker Compose Services

### Port Allocation

Existing HCLS AI Factory uses ports 8001-8510. Imaging agent uses 8520+:

| Service | Port | Protocol |
|---|---|---|
| Orthanc (DICOM) | 4242 | DIMSE |
| Orthanc (DICOMweb) | 8042 | HTTP REST |
| PostgreSQL | 5432 | TCP |
| NIM LLM | 8520 | HTTP (OpenAI-compatible) |
| Embedding Service | 8521 | HTTP REST |
| DICOM Listener | 8522 | HTTP webhook |
| FHIR Publisher | 8523 | HTTP REST |
| Agent API | 8524 | HTTP REST |
| Streamlit Portal | 8525 | HTTP |
| Prometheus | 9099 | HTTP |
| Grafana | 3000 | HTTP |
| DCGM Exporter | 9400 | HTTP |

### docker-compose.yml

```yaml
version: "3.8"

services:
  # ── Data Layer ────────────────────────────────────────────
  orthanc:
    image: orthancteam/orthanc:24.1.2
    container_name: imaging-orthanc
    ports:
      - "4242:4242"   # DIMSE
      - "8042:8042"   # DICOMweb REST
    volumes:
      - orthanc_data:/var/lib/orthanc/db
      - ./config/orthanc.json:/etc/orthanc/orthanc.json:ro
      - ./config/scripts:/etc/orthanc/scripts:ro
    environment:
      - ORTHANC__DICOM_AET=IMAGING_AGENT
      - ORTHANC__DICOM_PORT=4242
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8042/system"]
      interval: 30s
      timeout: 10s
      retries: 3

  postgres:
    image: pgvector/pgvector:pg16
    container_name: imaging-postgres
    ports:
      - "5432:5432"
    volumes:
      - postgres_data:/var/lib/postgresql/data
      - ./db/init.sql:/docker-entrypoint-initdb.d/init.sql:ro
    environment:
      - POSTGRES_USER=${POSTGRES_USER:-imaging}
      - POSTGRES_PASSWORD=${POSTGRES_PASSWORD:-imaging_secret}
      - POSTGRES_DB=${POSTGRES_DB:-imaging_agent}
    restart: unless-stopped
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U imaging"]
      interval: 10s
      timeout: 5s
      retries: 5

  # ── Execution Layer ───────────────────────────────────────
  nim-llm:
    image: nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark
    container_name: imaging-nim-llm
    ports:
      - "8520:8000"
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
    environment:
      - NGC_API_KEY=${NGC_API_KEY}
    volumes:
      - nim_cache:/opt/nim/.cache
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/v1/health/ready"]
      interval: 30s
      timeout: 10s
      retries: 5

  embedding-service:
    build:
      context: ./services/embedding_service
      dockerfile: Dockerfile
    container_name: imaging-embedding
    ports:
      - "8521:8000"
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
    environment:
      - POSTGRES_URL=postgresql://${POSTGRES_USER:-imaging}:${POSTGRES_PASSWORD:-imaging_secret}@postgres:5432/${POSTGRES_DB:-imaging_agent}
      - MODEL_NAME=microsoft/BiomedCLIP-PubMedBERT_256-vit_base_patch16_224
    depends_on:
      postgres:
        condition: service_healthy
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  dicom-listener:
    build:
      context: ./services/dicom_listener
      dockerfile: Dockerfile
    container_name: imaging-dicom-listener
    ports:
      - "8522:8000"
    environment:
      - ORTHANC_URL=http://orthanc:8042
      - POSTGRES_URL=postgresql://${POSTGRES_USER:-imaging}:${POSTGRES_PASSWORD:-imaging_secret}@postgres:5432/${POSTGRES_DB:-imaging_agent}
      - NEXTFLOW_API_URL=http://localhost:8080
    depends_on:
      orthanc:
        condition: service_healthy
      postgres:
        condition: service_healthy
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  fhir-publisher:
    build:
      context: ./services/fhir_publisher
      dockerfile: Dockerfile
    container_name: imaging-fhir-publisher
    ports:
      - "8523:8000"
    environment:
      - FHIR_SERVER_URL=${FHIR_SERVER_URL:-http://localhost:8080/fhir}
      - POSTGRES_URL=postgresql://${POSTGRES_USER:-imaging}:${POSTGRES_PASSWORD:-imaging_secret}@postgres:5432/${POSTGRES_DB:-imaging_agent}
    depends_on:
      postgres:
        condition: service_healthy
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  agent:
    build:
      context: .
      dockerfile: agent/Dockerfile
    container_name: imaging-agent
    ports:
      - "8524:8000"
    environment:
      - NIM_LLM_URL=http://nim-llm:8000/v1
      - POSTGRES_URL=postgresql://${POSTGRES_USER:-imaging}:${POSTGRES_PASSWORD:-imaging_secret}@postgres:5432/${POSTGRES_DB:-imaging_agent}
      - ORTHANC_URL=http://orthanc:8042
      - EMBEDDING_URL=http://embedding-service:8000
    depends_on:
      nim-llm:
        condition: service_healthy
      postgres:
        condition: service_healthy
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  portal:
    build:
      context: ./services/portal
      dockerfile: Dockerfile
    container_name: imaging-portal
    ports:
      - "8525:8501"
    environment:
      - AGENT_URL=http://agent:8000
      - ORTHANC_URL=http://orthanc:8042
      - POSTGRES_URL=postgresql://${POSTGRES_USER:-imaging}:${POSTGRES_PASSWORD:-imaging_secret}@postgres:5432/${POSTGRES_DB:-imaging_agent}
    depends_on:
      - agent
      - orthanc
      - postgres
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8501/healthz"]
      interval: 30s
      timeout: 10s
      retries: 3

  # ── Monitoring ────────────────────────────────────────────
  dcgm-exporter:
    image: nvcr.io/nvidia/k8s/dcgm-exporter:3.3.5-3.4.0-ubuntu22.04
    container_name: imaging-dcgm
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: all
              capabilities: [gpu]
    ports:
      - "9400:9400"
    restart: unless-stopped

  prometheus:
    image: prom/prometheus:v2.48.0
    container_name: imaging-prometheus
    ports:
      - "9099:9090"
    volumes:
      - ./config/prometheus.yml:/etc/prometheus/prometheus.yml:ro
      - prometheus_data:/prometheus
    command:
      - '--config.file=/etc/prometheus/prometheus.yml'
      - '--storage.tsdb.retention.time=30d'
    restart: unless-stopped
    depends_on:
      - dcgm-exporter

  grafana:
    image: grafana/grafana:10.2.2
    container_name: imaging-grafana
    ports:
      - "3000:3000"
    environment:
      - GF_SECURITY_ADMIN_USER=${GRAFANA_USER:-admin}
      - GF_SECURITY_ADMIN_PASSWORD=${GRAFANA_PASSWORD:-changeme}
    volumes:
      - grafana_data:/var/lib/grafana
      - ./config/grafana/dashboards:/var/lib/grafana/dashboards:ro
    restart: unless-stopped
    depends_on:
      - prometheus

volumes:
  orthanc_data:
  postgres_data:
  nim_cache:
  prometheus_data:
  grafana_data:

networks:
  default:
    name: imaging-agent-network
```

---

## 5. PostgreSQL + pgvector Schema

### init.sql

```sql
-- Enable pgvector extension
CREATE EXTENSION IF NOT EXISTS vector;

-- ── Studies ────────────────────────────────────────────────
CREATE TABLE studies (
    id                  SERIAL PRIMARY KEY,
    study_instance_uid  TEXT UNIQUE NOT NULL,
    patient_id          TEXT NOT NULL,
    patient_name        TEXT,
    study_date          DATE NOT NULL,
    study_description   TEXT,
    modality            TEXT NOT NULL,          -- CT, MR, CR, DX
    accession_number    TEXT,
    referring_physician TEXT,
    body_part           TEXT,                   -- HEAD, CHEST, BRAIN
    num_series          INT DEFAULT 0,
    num_instances       INT DEFAULT 0,
    orthanc_id          TEXT,                   -- Orthanc internal ID
    status              TEXT DEFAULT 'received', -- received, processing, completed, failed
    created_at          TIMESTAMPTZ DEFAULT NOW(),
    updated_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_studies_patient ON studies(patient_id);
CREATE INDEX idx_studies_modality ON studies(modality);
CREATE INDEX idx_studies_date ON studies(study_date);
CREATE INDEX idx_studies_status ON studies(status);

-- ── Series ─────────────────────────────────────────────────
CREATE TABLE series (
    id                  SERIAL PRIMARY KEY,
    series_instance_uid TEXT UNIQUE NOT NULL,
    study_id            INT REFERENCES studies(id) ON DELETE CASCADE,
    series_number       INT,
    series_description  TEXT,
    modality            TEXT NOT NULL,
    num_instances       INT DEFAULT 0,
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_series_study ON series(study_id);

-- ── Findings ───────────────────────────────────────────────
CREATE TABLE findings (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    workflow         TEXT NOT NULL,             -- ct_head_hemorrhage, ct_chest_nodule, cxr_findings, mri_ms_lesion
    finding_type    TEXT NOT NULL,             -- hemorrhage, nodule, pneumothorax, lesion, etc.
    finding_code    TEXT,                       -- SNOMED CT code
    location        TEXT,                       -- anatomical location
    laterality      TEXT,                       -- left, right, bilateral, midline
    severity        TEXT,                       -- critical, urgent, moderate, routine
    confidence      FLOAT NOT NULL CHECK (confidence >= 0 AND confidence <= 1),
    is_positive     BOOLEAN DEFAULT TRUE,
    details         JSONB DEFAULT '{}',        -- workflow-specific structured data
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_findings_study ON findings(study_id);
CREATE INDEX idx_findings_workflow ON findings(workflow);
CREATE INDEX idx_findings_severity ON findings(severity);
CREATE INDEX idx_findings_type ON findings(finding_type);
CREATE INDEX idx_findings_details ON findings USING GIN(details);

-- ── Measurements ───────────────────────────────────────────
CREATE TABLE measurements (
    id              SERIAL PRIMARY KEY,
    finding_id      INT REFERENCES findings(id) ON DELETE CASCADE,
    measurement_type TEXT NOT NULL,            -- volume, diameter, shift, count, ratio, doubling_time
    value           FLOAT NOT NULL,
    unit            TEXT NOT NULL,             -- mL, mm, mm3, days, ratio
    reference_range TEXT,                      -- e.g., "< 5 mm"
    flag            TEXT,                      -- normal, elevated, critical
    prior_value     FLOAT,                    -- from longitudinal comparison
    prior_date      DATE,
    delta_percent   FLOAT,                    -- percentage change from prior
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_measurements_finding ON measurements(finding_id);
CREATE INDEX idx_measurements_type ON measurements(measurement_type);

-- ── Embeddings ─────────────────────────────────────────────
CREATE TABLE embeddings (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    finding_id      INT REFERENCES findings(id) ON DELETE SET NULL,
    level           TEXT NOT NULL,             -- study, series, lesion
    model_name      TEXT NOT NULL,             -- embedding model identifier
    embedding       vector(384) NOT NULL,      -- 384-dim embedding vector
    metadata        JSONB DEFAULT '{}',
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

-- HNSW index for fast approximate nearest neighbor search
CREATE INDEX idx_embeddings_hnsw ON embeddings
    USING hnsw (embedding vector_cosine_ops)
    WITH (m = 16, ef_construction = 64);

CREATE INDEX idx_embeddings_study ON embeddings(study_id);
CREATE INDEX idx_embeddings_level ON embeddings(level);

-- ── Provenance ─────────────────────────────────────────────
CREATE TABLE provenance (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    workflow         TEXT NOT NULL,
    model_id        TEXT NOT NULL,             -- e.g., hemorrhage-triage-v2.1
    model_version   TEXT NOT NULL,             -- e.g., 2.1.0
    model_arch      TEXT,                      -- e.g., 3D U-Net (MONAI)
    inference_params JSONB DEFAULT '{}',       -- precision, seed, thresholds
    input_uids      TEXT[] DEFAULT '{}',       -- DICOM SOP Instance UIDs processed
    duration_ms     INT,                       -- inference duration
    gpu_memory_mb   INT,                       -- peak GPU memory
    status          TEXT DEFAULT 'completed',  -- completed, failed, timeout
    error_message   TEXT,
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_provenance_study ON provenance(study_id);
CREATE INDEX idx_provenance_workflow ON provenance(workflow);
CREATE INDEX idx_provenance_model ON provenance(model_id);

-- ── Worklist Entries ───────────────────────────────────────
CREATE TABLE worklist_entries (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    finding_id      INT REFERENCES findings(id) ON DELETE SET NULL,
    urgency         TEXT NOT NULL,             -- critical, urgent, moderate, routine
    priority        TEXT NOT NULL,             -- P1, P2, P3, P4
    notification    TEXT,                      -- alert text
    routing         TEXT,                      -- department/specialist
    acknowledged    BOOLEAN DEFAULT FALSE,
    acknowledged_by TEXT,
    acknowledged_at TIMESTAMPTZ,
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_worklist_urgency ON worklist_entries(urgency);
CREATE INDEX idx_worklist_ack ON worklist_entries(acknowledged);

-- ── Helper Views ───────────────────────────────────────────

-- Active worklist (unacknowledged, ordered by priority)
CREATE VIEW active_worklist AS
SELECT
    w.id, w.urgency, w.priority, w.notification, w.routing,
    s.patient_id, s.patient_name, s.modality, s.study_date,
    f.finding_type, f.confidence,
    w.created_at
FROM worklist_entries w
JOIN studies s ON w.study_id = s.id
LEFT JOIN findings f ON w.finding_id = f.id
WHERE w.acknowledged = FALSE
ORDER BY
    CASE w.priority
        WHEN 'P1' THEN 1 WHEN 'P2' THEN 2
        WHEN 'P3' THEN 3 ELSE 4
    END,
    w.created_at ASC;

-- Study summary with finding counts
CREATE VIEW study_summary AS
SELECT
    s.id, s.study_instance_uid, s.patient_id, s.modality,
    s.study_date, s.status,
    COUNT(DISTINCT f.id) AS finding_count,
    COUNT(DISTINCT f.id) FILTER (WHERE f.severity = 'critical') AS critical_count,
    COUNT(DISTINCT f.id) FILTER (WHERE f.severity = 'urgent') AS urgent_count,
    MAX(f.confidence) AS max_confidence,
    s.created_at
FROM studies s
LEFT JOIN findings f ON f.study_id = s.id
GROUP BY s.id;
```

### Example Hybrid Queries

```sql
-- All Lung-RADS 4A+ findings
SELECT f.*, m.value AS volume_mm3, m.delta_percent
FROM findings f
JOIN measurements m ON m.finding_id = f.id AND m.measurement_type = 'volume'
WHERE f.workflow = 'ct_chest_nodule'
  AND f.details->>'lung_rads' IN ('4A', '4B', '4X')
ORDER BY m.value DESC;

-- 10 most similar CT chest studies (vector search)
SELECT s.study_instance_uid, s.patient_id, s.study_date,
       e.embedding <=> $1::vector AS distance
FROM embeddings e
JOIN studies s ON e.study_id = s.id
WHERE e.level = 'study'
  AND s.modality = 'CT'
  AND s.body_part = 'CHEST'
ORDER BY e.embedding <=> $1::vector
LIMIT 10;

-- Growing nodules AND similar phenotype (hybrid)
WITH growing_nodules AS (
    SELECT f.study_id, f.id AS finding_id
    FROM findings f
    JOIN measurements m ON m.finding_id = f.id
    WHERE f.workflow = 'ct_chest_nodule'
      AND m.measurement_type = 'doubling_time'
      AND m.value < 400
)
SELECT s.patient_id, s.study_date,
       e.embedding <=> $1::vector AS phenotype_distance
FROM growing_nodules gn
JOIN studies s ON gn.study_id = s.id
JOIN embeddings e ON e.study_id = s.id AND e.level = 'study'
ORDER BY e.embedding <=> $1::vector
LIMIT 10;
```

---

## 6. Pydantic Data Models

### src/models.py

```python
"""Pydantic data models for the Imaging Intelligence Agent."""

from __future__ import annotations

from datetime import date, datetime
from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field


# ── Enums ───────────────────────────────────────────────────

class Modality(str, Enum):
    CT = "CT"
    MR = "MR"
    CR = "CR"
    DX = "DX"

class BodyPart(str, Enum):
    HEAD = "HEAD"
    CHEST = "CHEST"
    BRAIN = "BRAIN"
    ABDOMEN = "ABDOMEN"
    MSK = "MSK"

class Severity(str, Enum):
    CRITICAL = "critical"
    URGENT = "urgent"
    MODERATE = "moderate"
    ROUTINE = "routine"

class Priority(str, Enum):
    P1 = "P1"
    P2 = "P2"
    P3 = "P3"
    P4 = "P4"

class WorkflowType(str, Enum):
    CT_HEAD_HEMORRHAGE = "ct_head_hemorrhage"
    CT_CHEST_NODULE = "ct_chest_nodule"
    CXR_FINDINGS = "cxr_findings"
    MRI_MS_LESION = "mri_ms_lesion"

class StudyStatus(str, Enum):
    RECEIVED = "received"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"

class DiseaseActivity(str, Enum):
    STABLE = "stable"
    ACTIVE = "active"
    HIGHLY_ACTIVE = "highly_active"

class LungRADS(str, Enum):
    CAT_1 = "1"
    CAT_2 = "2"
    CAT_3 = "3"
    CAT_4A = "4A"
    CAT_4B = "4B"
    CAT_4X = "4X"


# ── Core Models ─────────────────────────────────────────────

class Study(BaseModel):
    """Represents a DICOM imaging study."""
    study_instance_uid: str
    patient_id: str
    patient_name: Optional[str] = None
    study_date: date
    study_description: Optional[str] = None
    modality: Modality
    accession_number: Optional[str] = None
    referring_physician: Optional[str] = None
    body_part: Optional[BodyPart] = None
    num_series: int = 0
    num_instances: int = 0
    orthanc_id: Optional[str] = None
    status: StudyStatus = StudyStatus.RECEIVED


class Finding(BaseModel):
    """A clinical finding produced by a workflow."""
    study_id: int
    workflow: WorkflowType
    finding_type: str
    finding_code: Optional[str] = None         # SNOMED CT
    location: Optional[str] = None
    laterality: Optional[str] = None
    severity: Severity
    confidence: float = Field(ge=0.0, le=1.0)
    is_positive: bool = True
    details: dict[str, Any] = Field(default_factory=dict)


class Measurement(BaseModel):
    """A quantitative measurement associated with a finding."""
    finding_id: int
    measurement_type: str                      # volume, diameter, shift, count, doubling_time
    value: float
    unit: str                                  # mL, mm, mm3, days
    reference_range: Optional[str] = None
    flag: Optional[str] = None                 # normal, elevated, critical
    prior_value: Optional[float] = None
    prior_date: Optional[date] = None
    delta_percent: Optional[float] = None


class EmbeddingRecord(BaseModel):
    """An embedding vector for similarity search."""
    study_id: int
    finding_id: Optional[int] = None
    level: str                                 # study, series, lesion
    model_name: str
    embedding: list[float]                     # 384-dim
    metadata: dict[str, Any] = Field(default_factory=dict)


class ProvenanceBundle(BaseModel):
    """Audit trail for an inference run."""
    study_id: int
    workflow: WorkflowType
    model_id: str
    model_version: str
    model_arch: Optional[str] = None
    inference_params: dict[str, Any] = Field(default_factory=dict)
    input_uids: list[str] = Field(default_factory=list)
    duration_ms: Optional[int] = None
    gpu_memory_mb: Optional[int] = None
    status: str = "completed"
    error_message: Optional[str] = None


class WorklistEntry(BaseModel):
    """A triage entry for the radiologist worklist."""
    study_id: int
    finding_id: Optional[int] = None
    urgency: Severity
    priority: Priority
    notification: Optional[str] = None
    routing: Optional[str] = None


# ── Workflow-Specific Models ────────────────────────────────

class HemorrhageResult(BaseModel):
    """CT Head hemorrhage triage result."""
    detected: bool
    hemorrhage_type: Optional[str] = None      # subdural, epidural, subarachnoid, intraparenchymal
    location: Optional[str] = None
    volume_ml: Optional[float] = None
    max_thickness_mm: Optional[float] = None
    midline_shift_mm: Optional[float] = None
    midline_shift_direction: Optional[str] = None
    urgency: Severity = Severity.ROUTINE
    confidence: float = Field(ge=0.0, le=1.0)


class NoduleResult(BaseModel):
    """Single lung nodule result."""
    nodule_id: str
    location: str
    nodule_type: str                           # solid, ground-glass, part-solid
    long_axis_mm: float
    short_axis_mm: float
    volume_mm3: float
    lung_rads: LungRADS
    malignancy_risk: Optional[float] = None
    confidence: float = Field(ge=0.0, le=1.0)
    # Longitudinal
    prior_volume_mm3: Optional[float] = None
    volume_change_percent: Optional[float] = None
    doubling_time_days: Optional[float] = None


class CXRFindingResult(BaseModel):
    """Single CXR finding result."""
    finding_name: str                          # pneumothorax, consolidation, etc.
    detected: bool
    confidence: float = Field(ge=0.0, le=1.0)
    gradcam_region: Optional[str] = None
    clinical_significance: Optional[str] = None


class MSLesionResult(BaseModel):
    """MRI Brain MS lesion tracking result."""
    total_lesion_count: int
    total_lesion_volume_ml: float
    new_lesion_count: int = 0
    enlarging_lesion_count: int = 0
    enhancing_lesion_count: int = 0
    disease_activity: DiseaseActivity
    prior_lesion_count: Optional[int] = None
    prior_volume_ml: Optional[float] = None
    volume_change_percent: Optional[float] = None
```

---

## 7. Orthanc DICOM Server Configuration

### config/orthanc.json

```json
{
    "Name": "ImagingAgent",
    "DicomAet": "IMAGING_AGENT",
    "DicomPort": 4242,
    "HttpPort": 8042,
    "StorageDirectory": "/var/lib/orthanc/db",
    "IndexDirectory": "/var/lib/orthanc/db",
    "StorageCompression": false,
    "RemoteAccessAllowed": true,
    "AuthenticationEnabled": false,

    "DicomWeb": {
        "Enable": true,
        "Root": "/dicom-web/",
        "EnableWado": true,
        "WadoRoot": "/wado",
        "StudiesMetadata": "Full"
    },

    "Lua": {
        "Enable": true,
        "Scripts": ["/etc/orthanc/scripts/on-stable-study.lua"]
    },

    "StableAge": 10,

    "DicomModalities": {
        "PACS": ["PACS_AET", "pacs-server", 4242]
    }
}
```

### Lua Script for study.complete Events

Save as `config/scripts/on-stable-study.lua`:

```lua
-- Fires when a study has been stable (no new instances) for StableAge seconds
function OnStableStudy(studyId, tags, metadata)
    -- Notify the DICOM listener service via HTTP webhook
    -- Port 8000 is the container-internal port; mapped to host 8522 in docker-compose
    local url = "http://dicom-listener:8000/webhook/study-complete"
    local body = '{"orthanc_id": "' .. studyId .. '"}'

    HttpPost(url, body, {["Content-Type"] = "application/json"})

    PrintToLog("Study stable, webhook sent: " .. studyId)
end
```

### DICOMweb Python Client Usage

```python
"""DICOM operations using dicomweb-client."""

from dicomweb_client import DICOMwebClient

ORTHANC_URL = "http://orthanc:8042/dicom-web"

client = DICOMwebClient(url=ORTHANC_URL)


def query_studies(patient_id: str | None = None, modality: str | None = None):
    """QIDO-RS: Query studies."""
    search_filters = {}
    if patient_id:
        search_filters["PatientID"] = patient_id
    if modality:
        search_filters["ModalitiesInStudy"] = modality
    return client.search_for_studies(search_filters=search_filters)


def retrieve_series(study_uid: str):
    """WADO-RS: Retrieve all series metadata for a study."""
    return client.retrieve_study_metadata(study_instance_uid=study_uid)


def retrieve_instances(study_uid: str, series_uid: str):
    """WADO-RS: Retrieve DICOM instances as pydicom datasets."""
    return client.retrieve_series(
        study_instance_uid=study_uid,
        series_instance_uid=series_uid,
    )


def store_instances(datasets: list):
    """STOW-RS: Store DICOM instances (SR, SEG, GSPS)."""
    client.store_instances(datasets=datasets)
```

---

## 8. MONAI Deploy MAP Pattern

### Application Class Structure

Each MONAI Deploy Application Package (MAP) follows this pattern:

```python
"""MONAI Deploy Application Package — [Workflow Name]."""

import logging

from monai.deploy.core import Application, resource
from monai.deploy.core import DataPath, ExecutionContext, InputContext, OutputContext


@resource(cpu=4, gpu=1, memory="16Gi")
class WorkflowApp(Application):
    """MONAI Deploy application for [workflow]."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._logger = logging.getLogger(self.__class__.__name__)

    def compose(self):
        """Define the operator pipeline."""
        preprocess_op = PreprocessOperator()
        inference_op = InferenceOperator()
        postprocess_op = PostprocessOperator()

        # Chain operators: input → preprocess → inference → postprocess → output
        self.add_flow(preprocess_op, inference_op)
        self.add_flow(inference_op, postprocess_op)


if __name__ == "__main__":
    WorkflowApp().run()
```

### Operator Pattern

```python
"""Operators for [workflow] MAP."""

import numpy as np
import torch
from monai.deploy.core import (
    ExecutionContext,
    InputContext,
    OutputContext,
    Operator,
    OperatorSpec,
)
from monai.transforms import (
    Compose,
    EnsureChannelFirstd,
    LoadImaged,
    Orientationd,
    ScaleIntensityRanged,
    Spacingd,
)


class PreprocessOperator(Operator):
    """Load and preprocess DICOM data."""

    def setup(self, spec: OperatorSpec):
        spec.input("dicom_path")
        spec.output("preprocessed")

    def compute(self, op_input: InputContext, op_output: OutputContext, context: ExecutionContext):
        dicom_path = op_input.get("dicom_path")

        transforms = Compose([
            LoadImaged(keys=["image"], reader="PydicomReader"),
            EnsureChannelFirstd(keys=["image"]),
            Orientationd(keys=["image"], axcodes="RAS"),
            Spacingd(keys=["image"], pixdim=(1.0, 1.0, 1.0), mode="bilinear"),
            ScaleIntensityRanged(keys=["image"], a_min=-1024, a_max=3071,
                                 b_min=0.0, b_max=1.0, clip=True),
        ])

        data = transforms({"image": str(dicom_path)})
        op_output.set(data, "preprocessed")


class InferenceOperator(Operator):
    """Run model inference."""

    def setup(self, spec: OperatorSpec):
        spec.input("preprocessed")
        spec.output("predictions")

    def compute(self, op_input: InputContext, op_output: OutputContext, context: ExecutionContext):
        data = op_input.get("preprocessed")
        image = data["image"]

        # Load model from mounted weights
        model_path = context.models.get("model")
        model = torch.jit.load(str(model_path), map_location="cuda")
        model.eval()

        with torch.no_grad():
            input_tensor = torch.from_numpy(image).unsqueeze(0).to("cuda")
            output = model(input_tensor)

        op_output.set({"predictions": output, "image": image}, "predictions")
```

### MAP Dockerfile Pattern

```dockerfile
# ARM64-compatible base for DGX Spark
FROM nvcr.io/nvidia/pytorch:24.01-py3

WORKDIR /app

# Install MONAI Deploy SDK
RUN pip install --no-cache-dir \
    monai-deploy-app-sdk>=0.6.0 \
    monai>=1.3.0 \
    pydicom>=2.4.0 \
    highdicom>=0.22.0 \
    numpy \
    scipy

COPY app.py operators.py requirements.txt ./

# Model weights mounted at runtime
ENV MONAI_DEPLOY_MODEL_PATH=/models

ENTRYPOINT ["python", "app.py"]
```

### I/O Conventions

| Path | Purpose |
|---|---|
| `/var/holoscan/input/` | DICOM input directory (mounted by orchestrator) |
| `/var/holoscan/output/` | Output directory (SR, SEG, measurements JSON) |
| `/models/` | Model weights directory (mounted volume) |

---

## 9. CT Head Hemorrhage Workflow

### Pipeline Stages

1. **Preprocessing:** Load NCCT head DICOM → reorient RAS → resample to 1mm isotropic → window (W:80 L:40 for blood) → normalize [0,1]
2. **Detection/Segmentation:** 3D U-Net binary segmentation (hemorrhage vs. normal)
3. **Volume Estimation:** Count positive voxels × voxel volume
4. **Midline Shift:** Compute brain centerline from falx cerebri, measure max lateral displacement
5. **Classification:** Map volume + shift to urgency (critical / urgent / routine)
6. **Output:** Finding + Measurements + WorklistEntry + DICOM SR

### MONAI Transforms Pipeline

```python
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd, Orientationd,
    Spacingd, ScaleIntensityRanged, CropForegroundd,
    EnsureTyped,
)

CT_HEAD_PREPROCESS = Compose([
    LoadImaged(keys=["image"], reader="PydicomReader"),
    EnsureChannelFirstd(keys=["image"]),
    Orientationd(keys=["image"], axcodes="RAS"),
    Spacingd(keys=["image"], pixdim=(1.0, 1.0, 1.0), mode="bilinear"),
    ScaleIntensityRanged(
        keys=["image"],
        a_min=0, a_max=80,       # CT brain window (blood)
        b_min=0.0, b_max=1.0,
        clip=True,
    ),
    CropForegroundd(keys=["image"], source_key="image", margin=10),
    EnsureTyped(keys=["image"], dtype="float32"),
])
```

### Model Architecture

```python
from monai.networks.nets import UNet

model = UNet(
    spatial_dims=3,
    in_channels=1,
    out_channels=2,          # background, hemorrhage
    channels=(16, 32, 64, 128, 256),
    strides=(2, 2, 2, 2),
    num_res_units=2,
    norm="batch",
)
```

### Midline Shift Measurement

```python
import numpy as np
from scipy.ndimage import center_of_mass

def measure_midline_shift(segmentation: np.ndarray, voxel_spacing: tuple[float, ...]) -> dict:
    """Measure midline shift from hemorrhage segmentation.

    Args:
        segmentation: Binary mask (H, W, D) where 1 = hemorrhage
        voxel_spacing: (sx, sy, sz) in mm

    Returns:
        dict with shift_mm, shift_direction
    """
    # Brain midline is at the center of the axial plane (left-right axis)
    axial_center = segmentation.shape[0] / 2.0  # assuming RAS orientation, axis 0 = L-R

    if segmentation.sum() == 0:
        return {"shift_mm": 0.0, "shift_direction": "none"}

    # Center of mass of hemorrhage in the L-R axis
    com = center_of_mass(segmentation)
    shift_voxels = com[0] - axial_center
    shift_mm = abs(shift_voxels * voxel_spacing[0])
    direction = "rightward" if shift_voxels > 0 else "leftward"

    return {"shift_mm": round(shift_mm, 1), "shift_direction": direction}
```

### Urgency Classification

```python
def classify_urgency(volume_ml: float, shift_mm: float, thickness_mm: float) -> tuple[str, str]:
    """Classify urgency based on Brain Trauma Foundation thresholds.

    Returns:
        (severity, priority) e.g. ("critical", "P1")
    """
    # BTF surgical thresholds: thickness > 10mm OR shift > 5mm
    if volume_ml > 30 or shift_mm > 5 or thickness_mm > 10:
        return ("critical", "P1")
    elif volume_ml > 5:
        return ("urgent", "P2")
    else:
        return ("routine", "P4")
```

### Performance Target

- **End-to-end:** < 90 seconds
- **Sensitivity:** > 95% for hemorrhage > 5 mL
- **Validated on:** RSNA ICH Dataset (752K slices)

---

## 10. CT Chest Lung Nodule Workflow

### Pipeline Stages

1. **Preprocessing:** Load CT chest DICOM → reorient RAS → resample to 1mm isotropic → lung window (W:1500 L:-600) → normalize
2. **Detection:** RetinaNet detects candidate nodules with bounding boxes
3. **Segmentation:** Per-nodule 3D U-Net (SegResNet) segmentation within each bounding box
4. **Volumetrics:** Voxel counting × voxel spacing product
5. **Prior Retrieval:** Query PostgreSQL for prior CT chest → retrieve from Orthanc
6. **Longitudinal Registration:** SyN diffeomorphic registration (ANTsPy) to align current and prior
7. **Volume Doubling Time:** `VDT = (Δt × ln2) / ln(V2/V1)`
8. **Lung-RADS Classification:** Rule-based assignment per ACR Lung-RADS v2022
9. **Malignancy Risk Scoring:** Composite score from VDT, morphology, location, patient factors
10. **Cross-modal Trigger:** If Lung-RADS 4B+ → trigger Parabricks genomics pipeline

### Detection Model

```python
from monai.networks.nets import RetinaNet
from monai.apps.detection.networks.retinanet_detector import RetinaNetDetector

detector = RetinaNetDetector(
    network=RetinaNet(
        spatial_dims=3,
        num_classes=1,         # nodule
        num_anchors=6,
        feature_extractor=resnet_fpn_feature_extractor(
            backbone="resnet50",
            spatial_dims=3,
            pretrained=False,
        ),
    ),
    anchor_generator=AnchorGeneratorWithAspectRatio(
        sizes=((4,), (6,), (8,), (12,), (16,)),
        aspect_ratios=((0.5, 1.0, 2.0),) * 5,
    ),
    score_thresh=0.5,
    nms_thresh=0.3,
)
```

### Segmentation Model

```python
from monai.networks.nets import SegResNet

seg_model = SegResNet(
    spatial_dims=3,
    in_channels=1,
    out_channels=2,           # background, nodule
    init_filters=16,
    blocks_down=(1, 2, 2, 4),
    blocks_up=(1, 1, 1),
    norm="batch",
)
```

### Volume Doubling Time

```python
import math
from datetime import date

def calculate_vdt(
    current_volume: float,
    prior_volume: float,
    current_date: date,
    prior_date: date,
) -> float | None:
    """Calculate volume doubling time in days.

    Formula: VDT = (Δt × ln2) / ln(V2/V1)
    """
    if prior_volume <= 0 or current_volume <= prior_volume:
        return None

    delta_days = (current_date - prior_date).days
    if delta_days <= 0:
        return None

    vdt = (delta_days * math.log(2)) / math.log(current_volume / prior_volume)
    return round(vdt, 1)
```

### Lung-RADS Classification

```python
def assign_lung_rads(
    nodule_type: str,
    long_axis_mm: float,
    volume_mm3: float,
    is_new: bool = False,
    prior_lung_rads: str | None = None,
    vdt_days: float | None = None,
) -> str:
    """Assign Lung-RADS category per ACR Lung-RADS v2022.

    Returns one of: '1', '2', '3', '4A', '4B', '4X'
    """
    # Solid nodules
    if nodule_type == "solid":
        if long_axis_mm < 4:
            return "1"
        elif long_axis_mm < 6:
            return "2"
        elif long_axis_mm < 8:
            return "3"
        elif long_axis_mm < 15:
            category = "4A"
        else:
            category = "4B"

    # Ground-glass nodules
    elif nodule_type == "ground-glass":
        if long_axis_mm < 6:
            return "1"
        elif long_axis_mm < 20:
            return "2"
        elif long_axis_mm < 30:
            return "3"
        else:
            category = "4A"

    # Part-solid nodules
    elif nodule_type == "part-solid":
        if long_axis_mm < 6:
            return "2"
        elif long_axis_mm < 8:
            return "3"
        else:
            category = "4A"

    else:
        category = "3"

    # Upgrade for growth
    if vdt_days is not None and vdt_days < 400:
        if category in ("2", "3"):
            category = "4A"
        elif category == "4A":
            category = "4B"

    # 4X: additional features (spiculation, lymphadenopathy, etc.)
    # Handled separately in the full implementation

    return category
```

### Genomics Pipeline Trigger

```python
def check_genomics_trigger(lung_rads: str, study_id: int) -> bool:
    """Trigger Parabricks genomics pipeline if Lung-RADS 4B+."""
    if lung_rads in ("4B", "4X"):
        # POST to Nextflow API or message queue
        trigger_genomics_pipeline(study_id=study_id, reason=f"Lung-RADS {lung_rads}")
        return True
    return False
```

---

## 11. CXR Rapid Findings Workflow

### Pipeline Stages

1. **Preprocessing:** Load CXR DICOM → resize to 224×224 → normalize (ImageNet stats)
2. **Multi-label Classification:** DenseNet-121 with 5 output heads
3. **GradCAM Heatmap:** Generate localization map for each positive finding
4. **Thresholding:** Per-class confidence thresholds
5. **Output:** Findings + GradCAM images (Secondary Capture) + GSPS overlays

### Classification Model

```python
from monai.networks.nets import DenseNet121

model = DenseNet121(
    spatial_dims=2,
    in_channels=1,         # grayscale CXR
    out_channels=5,        # pneumothorax, consolidation, pleural_effusion, cardiomegaly, fracture
)

FINDING_LABELS = [
    "pneumothorax",
    "consolidation",
    "pleural_effusion",
    "cardiomegaly",
    "fracture",
]

CONFIDENCE_THRESHOLDS = {
    "pneumothorax": 0.50,       # Lower threshold — high-risk finding
    "consolidation": 0.60,
    "pleural_effusion": 0.55,
    "cardiomegaly": 0.60,
    "fracture": 0.55,
}
```

### Preprocessing

```python
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd,
    Resized, ScaleIntensityd, EnsureTyped,
)

CXR_PREPROCESS = Compose([
    LoadImaged(keys=["image"], reader="PydicomReader"),
    EnsureChannelFirstd(keys=["image"]),
    Resized(keys=["image"], spatial_size=(224, 224), mode="bilinear"),
    ScaleIntensityd(keys=["image"]),  # normalize to [0, 1]
    EnsureTyped(keys=["image"], dtype="float32"),
])
```

### GradCAM Heatmap Generation

```python
import torch
import numpy as np
from monai.visualize import GradCAM

def generate_gradcam(
    model: torch.nn.Module,
    input_tensor: torch.Tensor,
    target_class: int,
) -> np.ndarray:
    """Generate GradCAM heatmap for a specific class.

    Args:
        model: DenseNet-121 model
        input_tensor: (1, 1, 224, 224) input
        target_class: class index (0-4)

    Returns:
        heatmap: (224, 224) numpy array, values in [0, 1]
    """
    cam = GradCAM(
        nn_module=model,
        target_layers="class_layers.relu",   # DenseNet-121 final ReLU
    )

    # GradCAM returns (B, 1, H, W) for 2D
    heatmap = cam(x=input_tensor, class_idx=target_class)
    return heatmap.squeeze().cpu().numpy()
```

### Multi-label Output Parsing

```python
import torch

def parse_cxr_predictions(
    logits: torch.Tensor,
    thresholds: dict[str, float] = CONFIDENCE_THRESHOLDS,
) -> list[CXRFindingResult]:
    """Parse multi-label model output into findings."""
    probs = torch.sigmoid(logits).squeeze().cpu().numpy()
    results = []

    for i, label in enumerate(FINDING_LABELS):
        detected = probs[i] >= thresholds[label]
        results.append(CXRFindingResult(
            finding_name=label,
            detected=detected,
            confidence=round(float(probs[i]), 2),
        ))

    return results
```

### Performance Target

- **End-to-end:** < 30 seconds
- **Pneumothorax sensitivity:** > 95%
- **Validated on:** CheXpert + MIMIC-CXR (601K images)

---

## 12. MRI Brain MS Lesion Workflow

### Pipeline Stages

1. **Preprocessing:** Load MRI FLAIR DICOM → reorient RAS → resample to 1mm isotropic → normalize (z-score)
2. **Segmentation:** 3D U-Net on FLAIR for white matter lesion segmentation
3. **Connected Components:** Label individual lesions, count them
4. **Volumetrics:** Per-lesion and total volume measurement
5. **Prior Retrieval:** Query PostgreSQL for prior MRI brain → retrieve from Orthanc
6. **Registration:** Affine + SyN registration to align current and prior FLAIR
7. **Lesion Matching:** Overlap analysis to identify new vs. stable vs. enlarging lesions
8. **Disease Activity:** Classify as stable / active / highly active

### Segmentation Model

```python
from monai.networks.nets import UNet

ms_lesion_model = UNet(
    spatial_dims=3,
    in_channels=1,              # FLAIR
    out_channels=2,             # background, lesion
    channels=(32, 64, 128, 256),
    strides=(2, 2, 2),
    num_res_units=2,
    norm="batch",
)
```

### Preprocessing

```python
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd, Orientationd,
    Spacingd, NormalizeIntensityd, CropForegroundd, EnsureTyped,
)

MRI_MS_PREPROCESS = Compose([
    LoadImaged(keys=["image"], reader="PydicomReader"),
    EnsureChannelFirstd(keys=["image"]),
    Orientationd(keys=["image"], axcodes="RAS"),
    Spacingd(keys=["image"], pixdim=(1.0, 1.0, 1.0), mode="bilinear"),
    NormalizeIntensityd(keys=["image"], nonzero=True, channel_wise=True),
    CropForegroundd(keys=["image"], source_key="image", margin=10),
    EnsureTyped(keys=["image"], dtype="float32"),
])
```

### Connected Component Analysis

```python
import numpy as np
from scipy.ndimage import label as scipy_label

def analyze_lesions(
    segmentation: np.ndarray,
    voxel_spacing: tuple[float, float, float],
) -> list[dict]:
    """Identify individual lesions and measure each.

    Args:
        segmentation: Binary mask (H, W, D) where 1 = lesion
        voxel_spacing: (sx, sy, sz) in mm

    Returns:
        List of lesion dicts with id, volume_mm3, centroid
    """
    labeled_array, num_features = scipy_label(segmentation)
    voxel_volume_mm3 = float(np.prod(voxel_spacing))

    lesions = []
    for i in range(1, num_features + 1):
        mask = labeled_array == i
        volume_mm3 = float(mask.sum()) * voxel_volume_mm3
        centroid = np.array(np.where(mask)).mean(axis=1)

        lesions.append({
            "lesion_id": f"L{i}",
            "volume_mm3": round(volume_mm3, 1),
            "centroid_voxel": centroid.tolist(),
            "centroid_mm": (centroid * np.array(voxel_spacing)).tolist(),
        })

    return lesions
```

### Spatial Registration (ANTsPy)

```python
import ants

def register_mri(
    fixed_path: str,
    moving_path: str,
) -> tuple:
    """Register prior MRI to current using SyN diffeomorphic registration.

    Returns:
        (warped_image, transform) — transform can be applied to prior lesion masks
    """
    fixed = ants.image_read(fixed_path)
    moving = ants.image_read(moving_path)

    result = ants.registration(
        fixed=fixed,
        moving=moving,
        type_of_transform="SyNRA",   # Rigid + Affine + SyN
        verbose=False,
    )

    return result["warpedmovout"], result["fwdtransforms"]
```

### Lesion Matching and Disease Activity

```python
def match_lesions(
    current_lesions: list[dict],
    prior_lesions_warped: list[dict],
    overlap_threshold: float = 0.3,
) -> dict:
    """Match current lesions to prior (after registration).

    Returns:
        {"new": [...], "stable": [...], "enlarging": [...]}
    """
    matched = {"new": [], "stable": [], "enlarging": []}
    prior_matched = set()

    for curr in current_lesions:
        best_overlap = 0.0
        best_prior = None
        for j, prior in enumerate(prior_lesions_warped):
            overlap = compute_dice(curr["mask"], prior["mask"])
            if overlap > best_overlap:
                best_overlap = overlap
                best_prior = j

        if best_overlap >= overlap_threshold and best_prior is not None:
            prior_matched.add(best_prior)
            vol_change = (curr["volume_mm3"] - prior_lesions_warped[best_prior]["volume_mm3"])
            if vol_change > 50:  # > 50 mm3 growth threshold
                matched["enlarging"].append(curr)
            else:
                matched["stable"].append(curr)
        else:
            matched["new"].append(curr)

    return matched


def classify_disease_activity(
    new_count: int,
    enlarging_count: int,
    enhancing_count: int,
) -> str:
    """Classify MS disease activity.

    Returns: 'stable', 'active', or 'highly_active'
    """
    total_active = new_count + enlarging_count
    if total_active == 0 and enhancing_count == 0:
        return "stable"
    elif total_active >= 3 or enhancing_count >= 2:
        return "highly_active"
    else:
        return "active"
```

### Performance Target

- **End-to-end:** < 5 minutes (multi-stage with registration)
- **Validated on:** ISBI MS Challenge + institutional data (1,200 MRIs)

---

## 13. DICOM SR Output (highdicom)

### Creating a TID 1500 Measurement Report

```python
"""DICOM SR generation using highdicom."""

import highdicom as hd
from highdicom.sr import (
    CodedConcept,
    Comprehensive3DSR,
    FindingSite,
    MeasurementReport,
    Measurement,
    NumericMeasurement,
    QualitativeFinding,
    TrackingIdentifier,
)
from pydicom.sr.codedict import codes
import pydicom


def create_hemorrhage_sr(
    source_dataset: pydicom.Dataset,
    volume_ml: float,
    shift_mm: float,
    location: str,
    urgency: str,
    confidence: float,
) -> Comprehensive3DSR:
    """Create DICOM SR for CT Head hemorrhage finding."""

    # Tracking identifier for the finding
    tracking = TrackingIdentifier(
        identifier="hemorrhage-finding-1",
        uid=hd.UID(),
    )

    # Finding site
    finding_site = FindingSite(
        anatomic_location=CodedConcept(
            value="12738006",
            meaning="Brain",
            scheme_designator="SCT",
        ),
        laterality=CodedConcept(
            value="7771000",
            meaning="Left",
            scheme_designator="SCT",
        ),
    )

    # Measurements
    volume_measurement = NumericMeasurement(
        name=CodedConcept(
            value="118565006",
            meaning="Volume",
            scheme_designator="SCT",
        ),
        value=volume_ml,
        unit=CodedConcept(
            value="mL",
            meaning="milliliter",
            scheme_designator="UCUM",
        ),
    )

    shift_measurement = NumericMeasurement(
        name=CodedConcept(
            value="LP267085-3",
            meaning="Midline shift",
            scheme_designator="LN",
        ),
        value=shift_mm,
        unit=CodedConcept(
            value="mm",
            meaning="millimeter",
            scheme_designator="UCUM",
        ),
    )

    # Build the measurement report
    report = MeasurementReport(
        observation_context=hd.sr.ObservationContext(
            observer_person_context=None,
            observer_device_context=hd.sr.ObserverContext(
                observer_type=CodedConcept(
                    value="121007",
                    meaning="Device",
                    scheme_designator="DCM",
                ),
                observer_identifying_attributes=hd.sr.DeviceObserverIdentifyingAttributes(
                    uid=hd.UID(),
                    name="ImagingIntelligenceAgent",
                ),
            ),
        ),
        procedure_reported=CodedConcept(
            value="77477000",
            meaning="CT of head",
            scheme_designator="SCT",
        ),
        imaging_measurements=[volume_measurement, shift_measurement],
        title=CodedConcept(
            value="126000",
            meaning="Imaging Measurement Report",
            scheme_designator="DCM",
        ),
    )

    # Create Comprehensive 3D SR
    sr = Comprehensive3DSR(
        evidence=[source_dataset],
        content=report,
        series_number=100,
        series_instance_uid=hd.UID(),
        sop_instance_uid=hd.UID(),
        instance_number=1,
        manufacturer="HCLS AI Factory",
        institution_name="Imaging Intelligence Agent",
    )

    return sr
```

### Pushing SR to Orthanc

```python
from dicomweb_client import DICOMwebClient

def push_sr_to_pacs(sr_dataset, orthanc_url: str = "http://orthanc:8042/dicom-web"):
    """Store DICOM SR to Orthanc via STOW-RS."""
    client = DICOMwebClient(url=orthanc_url)
    client.store_instances(datasets=[sr_dataset])
```

---

## 14. FHIR DiagnosticReport Output

### Creating a FHIR DiagnosticReport

```python
"""FHIR DiagnosticReport generation."""

from datetime import datetime
from fhir.resources.diagnosticreport import DiagnosticReport
from fhir.resources.observation import Observation
from fhir.resources.codeableconcept import CodeableConcept
from fhir.resources.coding import Coding
from fhir.resources.reference import Reference
from fhir.resources.quantity import Quantity
from fhir.resources.bundle import Bundle, BundleEntry, BundleEntryRequest


def create_diagnostic_report(
    patient_id: str,
    study_uid: str,
    findings: list[dict],
    conclusion: str,
) -> Bundle:
    """Create a FHIR Bundle with DiagnosticReport and Observations."""

    observations = []
    observation_refs = []

    for i, finding in enumerate(findings):
        obs = Observation(
            status="final",
            code=CodeableConcept(
                coding=[Coding(
                    system="http://snomed.info/sct",
                    code=finding.get("snomed_code", "404684003"),
                    display=finding["finding_type"],
                )],
            ),
            valueString=finding.get("description", ""),
            interpretation=[CodeableConcept(
                coding=[Coding(
                    system="http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                    code="A" if finding.get("severity") in ("critical", "urgent") else "N",
                    display="Abnormal" if finding.get("severity") in ("critical", "urgent") else "Normal",
                )],
            )],
            subject=Reference(reference=f"Patient/{patient_id}"),
        )

        # Add quantitative components
        if "measurements" in finding:
            obs.component = []
            for meas in finding["measurements"]:
                obs.component.append({
                    "code": CodeableConcept(
                        coding=[Coding(
                            system="http://loinc.org",
                            code=meas.get("loinc_code", "LP6811-4"),
                            display=meas["type"],
                        )],
                    ),
                    "valueQuantity": Quantity(
                        value=meas["value"],
                        unit=meas["unit"],
                        system="http://unitsofmeasure.org",
                        code=meas["unit"],
                    ),
                })

        observations.append(obs)
        observation_refs.append(Reference(reference=f"Observation/{i}"))

    # DiagnosticReport
    report = DiagnosticReport(
        status="final",
        code=CodeableConcept(
            coding=[Coding(
                system="http://loinc.org",
                code="18748-4",
                display="Diagnostic imaging study",
            )],
        ),
        subject=Reference(reference=f"Patient/{patient_id}"),
        effectiveDateTime=datetime.now().isoformat(),
        conclusion=conclusion,
        result=observation_refs,
        imagingStudy=[Reference(
            reference=f"ImagingStudy/{study_uid}",
        )],
    )

    # Bundle everything
    entries = []
    entries.append(BundleEntry(
        resource=report,
        request=BundleEntryRequest(method="POST", url="DiagnosticReport"),
    ))
    for obs in observations:
        entries.append(BundleEntry(
            resource=obs,
            request=BundleEntryRequest(method="POST", url="Observation"),
        ))

    return Bundle(type="transaction", entry=entries)
```

---

## 15. LangGraph Agent Architecture

### Agent State

```python
"""Agent state definition."""

from typing import Annotated, TypedDict

from langgraph.graph.message import add_messages


class AgentState(TypedDict):
    """State shared across all agent nodes."""
    messages: Annotated[list, add_messages]
    study_id: int
    study_uid: str
    patient_id: str
    modality: str
    findings: list[dict]
    measurements: list[dict]
    worklist_entries: list[dict]
    prior_studies: list[dict]
    similar_studies: list[dict]
    rag_context: str
    report_text: str
    triage_complete: bool
    longitudinal_complete: bool
```

### Graph Definition

```python
"""LangGraph clinical reasoning agent."""

from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver

from agent.state import AgentState
from agent.nodes import (
    triage_node,
    longitudinal_node,
    population_node,
    report_node,
)


def build_agent_graph() -> StateGraph:
    """Build the clinical reasoning agent graph."""

    graph = StateGraph(AgentState)

    # Add nodes
    graph.add_node("triage", triage_node)
    graph.add_node("longitudinal", longitudinal_node)
    graph.add_node("population", population_node)
    graph.add_node("report", report_node)

    # Edges
    graph.add_edge(START, "triage")
    graph.add_conditional_edges(
        "triage",
        route_after_triage,
        {
            "longitudinal": "longitudinal",
            "report": "report",
        },
    )
    graph.add_edge("longitudinal", "population")
    graph.add_edge("population", "report")
    graph.add_edge("report", END)

    return graph.compile(checkpointer=MemorySaver())


def route_after_triage(state: AgentState) -> str:
    """Route based on finding severity — critical/urgent go to longitudinal analysis."""
    severities = [f.get("severity") for f in state.get("findings", [])]
    if "critical" in severities or "urgent" in severities:
        return "longitudinal"
    return "report"
```

### Agent Nodes

```python
"""Agent node implementations."""

from langchain_openai import ChatOpenAI
from agent.state import AgentState
from agent.tools import query_findings, search_similar_studies, get_prior_measurements


# LLM client pointing to local NIM
llm = ChatOpenAI(
    base_url="http://nim-llm:8000/v1",
    api_key="not-needed",          # NIM uses NGC_API_KEY from env
    model="meta-llama3-8b-instruct",
    temperature=0.1,
)


def triage_node(state: AgentState) -> dict:
    """Classify findings by urgency and create worklist entries."""
    findings = state["findings"]
    worklist = []

    for finding in findings:
        severity = finding.get("severity", "routine")
        priority_map = {
            "critical": "P1",
            "urgent": "P2",
            "moderate": "P3",
            "routine": "P4",
        }
        worklist.append({
            "study_id": state["study_id"],
            "finding_id": finding.get("id"),
            "urgency": severity,
            "priority": priority_map.get(severity, "P4"),
            "notification": f"Alert: {finding['finding_type']} — {severity}",
            "routing": determine_routing(finding),
        })

    return {"worklist_entries": worklist, "triage_complete": True}


def longitudinal_node(state: AgentState) -> dict:
    """Retrieve prior studies and compute longitudinal deltas."""
    prior_studies = get_prior_measurements.invoke({
        "patient_id": state["patient_id"],
        "modality": state["modality"],
    })
    return {"prior_studies": prior_studies, "longitudinal_complete": True}


def population_node(state: AgentState) -> dict:
    """Find similar studies via embedding search."""
    similar = search_similar_studies.invoke({
        "study_id": state["study_id"],
        "limit": 10,
    })
    return {"similar_studies": similar}


def report_node(state: AgentState) -> dict:
    """Generate evidence-grounded clinical summary via RAG."""
    # Build context from findings, priors, similar cases
    context_parts = []
    for f in state["findings"]:
        context_parts.append(f"Finding: {f['finding_type']} — {f.get('severity')} (confidence: {f.get('confidence')})")
    for p in state.get("prior_studies", []):
        context_parts.append(f"Prior: {p}")

    context = "\n".join(context_parts)

    response = llm.invoke(
        f"You are a radiology AI assistant. Summarize the following imaging findings "
        f"into a structured clinical report. Ground your recommendations in ACR guidelines.\n\n"
        f"{context}"
    )

    return {"report_text": response.content, "rag_context": context}


def determine_routing(finding: dict) -> str:
    """Determine specialist routing based on finding type."""
    routing_map = {
        "hemorrhage": "Neurosurgery",
        "nodule": "Pulmonology / Interventional Radiology",
        "pneumothorax": "Emergency Medicine / Thoracic Surgery",
        "ms_lesion": "Neurology / MS Clinic",
    }
    return routing_map.get(finding.get("finding_type", ""), "Radiology")
```

### MCP Tool Definitions

```python
"""MCP tool definitions for agent."""

from langchain_core.tools import tool
import psycopg2
from pgvector.psycopg2 import register_vector

from src.config import POSTGRES_URL


@tool
def query_findings(workflow: str, severity: str | None = None, limit: int = 50) -> list[dict]:
    """Query structured findings from the imaging database.

    Args:
        workflow: Workflow type (ct_head_hemorrhage, ct_chest_nodule, cxr_findings, mri_ms_lesion)
        severity: Optional severity filter (critical, urgent, moderate, routine)
        limit: Maximum results to return
    """
    conn = psycopg2.connect(POSTGRES_URL)
    register_vector(conn)
    cur = conn.cursor()

    query = "SELECT * FROM findings WHERE workflow = %s"
    params = [workflow]
    if severity:
        query += " AND severity = %s"
        params.append(severity)
    query += " ORDER BY created_at DESC LIMIT %s"
    params.append(limit)

    cur.execute(query, params)
    columns = [desc[0] for desc in cur.description]
    results = [dict(zip(columns, row)) for row in cur.fetchall()]
    conn.close()
    return results


@tool
def search_similar_studies(study_id: int, limit: int = 10) -> list[dict]:
    """Find similar studies using embedding vector search.

    Args:
        study_id: Source study ID
        limit: Number of similar studies to return
    """
    conn = psycopg2.connect(POSTGRES_URL)
    register_vector(conn)
    cur = conn.cursor()

    # Get the source embedding
    cur.execute(
        "SELECT embedding FROM embeddings WHERE study_id = %s AND level = 'study' LIMIT 1",
        [study_id],
    )
    row = cur.fetchone()
    if not row:
        conn.close()
        return []

    source_embedding = row[0]

    # Find nearest neighbors
    cur.execute(
        """SELECT s.study_instance_uid, s.patient_id, s.study_date, s.modality,
                  e.embedding <=> %s::vector AS distance
           FROM embeddings e
           JOIN studies s ON e.study_id = s.id
           WHERE e.level = 'study' AND e.study_id != %s
           ORDER BY e.embedding <=> %s::vector
           LIMIT %s""",
        [source_embedding, study_id, source_embedding, limit],
    )
    columns = [desc[0] for desc in cur.description]
    results = [dict(zip(columns, row)) for row in cur.fetchall()]
    conn.close()
    return results


@tool
def get_prior_measurements(patient_id: str, modality: str) -> list[dict]:
    """Retrieve prior study measurements for longitudinal comparison.

    Args:
        patient_id: Patient identifier
        modality: Imaging modality (CT, MR, CR)
    """
    conn = psycopg2.connect(POSTGRES_URL)
    cur = conn.cursor()

    cur.execute(
        """SELECT s.study_date, f.finding_type, f.location,
                  m.measurement_type, m.value, m.unit
           FROM studies s
           JOIN findings f ON f.study_id = s.id
           JOIN measurements m ON m.finding_id = f.id
           WHERE s.patient_id = %s AND s.modality = %s
           ORDER BY s.study_date DESC""",
        [patient_id, modality],
    )
    columns = [desc[0] for desc in cur.description]
    results = [dict(zip(columns, row)) for row in cur.fetchall()]
    conn.close()
    return results
```

---

## 16. NIM LLM Integration

### Deployment on DGX Spark

NIM on DGX Spark uses ARM64 containers with the `-dgx-spark` tag suffix:

```yaml
# In docker-compose.yml
nim-llm:
  image: nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark
  environment:
    - NGC_API_KEY=${NGC_API_KEY}
  deploy:
    resources:
      reservations:
        devices:
          - driver: nvidia
            count: 1
            capabilities: [gpu]
```

### OpenAI-Compatible API Usage

```python
"""NIM LLM client — OpenAI-compatible API."""

from openai import OpenAI

nim_client = OpenAI(
    base_url="http://nim-llm:8000/v1",
    api_key="not-needed",
)


def generate_clinical_summary(
    findings_context: str,
    guidelines_context: str,
    max_tokens: int = 1024,
) -> str:
    """Generate evidence-grounded clinical summary via RAG."""
    response = nim_client.chat.completions.create(
        model="meta-llama3-8b-instruct",
        messages=[
            {
                "role": "system",
                "content": (
                    "You are a radiology AI assistant. Provide structured clinical "
                    "summaries grounded in imaging evidence and ACR guidelines. "
                    "Always cite specific measurements and confidence scores."
                ),
            },
            {
                "role": "user",
                "content": (
                    f"## Imaging Findings\n{findings_context}\n\n"
                    f"## Relevant Guidelines\n{guidelines_context}\n\n"
                    "Generate a structured clinical summary with recommendations."
                ),
            },
        ],
        temperature=0.1,
        max_tokens=max_tokens,
    )
    return response.choices[0].message.content
```

### RAG Pipeline Pattern

```python
def rag_pipeline(study_id: int, findings: list[dict]) -> str:
    """Full RAG pipeline: retrieve evidence → construct prompt → generate."""

    # 1. Retrieve relevant guidelines from pgvector
    guidelines = search_guidelines(findings)

    # 2. Retrieve prior measurements for longitudinal context
    priors = get_prior_measurements(study_id)

    # 3. Retrieve similar cases for outcome context
    similar = search_similar_studies(study_id, limit=5)

    # 4. Construct context
    context = format_rag_context(findings, guidelines, priors, similar)

    # 5. Generate grounded response
    summary = generate_clinical_summary(
        findings_context=format_findings(findings),
        guidelines_context=context,
    )

    return summary
```

---

## 17. Embedding Service

### Embedding Generation

```python
"""Study/series/lesion embedding service."""

import numpy as np
import torch
from fastapi import FastAPI
from transformers import AutoModel, AutoTokenizer
from pydantic import BaseModel

app = FastAPI(title="Imaging Embedding Service")

# Load medical vision encoder
MODEL_NAME = "microsoft/BiomedCLIP-PubMedBERT_256-vit_base_patch16_224"
model = AutoModel.from_pretrained(MODEL_NAME).to("cuda").eval()


class EmbeddingRequest(BaseModel):
    study_id: int
    level: str          # study, series, lesion
    image_path: str     # path to DICOM or NIfTI


class EmbeddingResponse(BaseModel):
    study_id: int
    level: str
    embedding: list[float]


@app.post("/embed", response_model=EmbeddingResponse)
async def generate_embedding(request: EmbeddingRequest):
    """Generate 384-dim embedding for an imaging study/series/lesion."""
    # Load and preprocess image
    image_tensor = load_and_preprocess(request.image_path)

    with torch.no_grad():
        embedding = model.get_image_features(image_tensor.to("cuda"))
        embedding = embedding.cpu().numpy().flatten()

    # Normalize to unit vector
    embedding = embedding / np.linalg.norm(embedding)

    return EmbeddingResponse(
        study_id=request.study_id,
        level=request.level,
        embedding=embedding.tolist(),
    )
```

### Storing in pgvector

```python
import psycopg2
from pgvector.psycopg2 import register_vector
import numpy as np

def store_embedding(study_id: int, level: str, embedding: list[float], model_name: str):
    """Store embedding in pgvector."""
    conn = psycopg2.connect(POSTGRES_URL)
    register_vector(conn)
    cur = conn.cursor()

    cur.execute(
        "INSERT INTO embeddings (study_id, level, model_name, embedding) VALUES (%s, %s, %s, %s)",
        [study_id, level, model_name, np.array(embedding)],
    )
    conn.commit()
    conn.close()
```

---

## 18. Nextflow Pipeline Orchestration

### main.nf

```groovy
nextflow.enable.dsl = 2

// Pipeline metadata
def pipelineVersion = '1.0.0'
def pipelineName = 'Imaging-Intelligence-Agent'

// Include workflow modules
include { CT_HEAD_HEMORRHAGE } from './modules/ct_head_hemorrhage'
include { CT_CHEST_LUNG_NODULE } from './modules/ct_chest_lung_nodule'
include { CXR_RAPID_FINDINGS } from './modules/cxr_rapid_findings'
include { MRI_BRAIN_MS_LESION } from './modules/mri_brain_ms_lesion'

// Input validation
def validateInputs() {
    if (!params.study_uid && !params.input_dir) {
        log.info "No study_uid or input_dir provided — running in demo mode"
    }
}

// Main workflow
workflow {
    validateInputs()

    if (params.workflow == 'ct_head') {
        Channel.of(params.study_uid).set { ch_study }
        CT_HEAD_HEMORRHAGE(ch_study)

    } else if (params.workflow == 'ct_chest') {
        Channel.of(params.study_uid).set { ch_study }
        CT_CHEST_LUNG_NODULE(ch_study)

    } else if (params.workflow == 'cxr') {
        Channel.of(params.study_uid).set { ch_study }
        CXR_RAPID_FINDINGS(ch_study)

    } else if (params.workflow == 'mri_brain') {
        Channel.of(params.study_uid).set { ch_study }
        MRI_BRAIN_MS_LESION(ch_study)

    } else {
        log.info "Unknown workflow: ${params.workflow}"
        log.info "Available: ct_head, ct_chest, cxr, mri_brain"
    }
}

// Completion handler
workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    log.info "Pipeline ${pipelineName} v${pipelineVersion}: ${status}"
    log.info "Duration: ${workflow.duration}"
}
```

### Workflow Module Example (CT Head)

```groovy
// modules/ct_head_hemorrhage.nf

process PREPROCESS {
    tag "$study_uid"
    label 'gpu'
    container 'imaging-agent/ct-head-hemorrhage:latest'

    input:
    val study_uid

    output:
    tuple val(study_uid), path("preprocessed/"), emit: preprocessed
    path "versions.yml", emit: versions

    script:
    """
    python -m maps.ct_head_hemorrhage.preprocess \
        --study-uid ${study_uid} \
        --orthanc-url ${params.orthanc_url} \
        --output-dir preprocessed/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        monai: \$(python -c "import monai; print(monai.__version__)")
    END_VERSIONS
    """
}

process INFERENCE {
    tag "$study_uid"
    label 'gpu'
    container 'imaging-agent/ct-head-hemorrhage:latest'

    input:
    tuple val(study_uid), path(preprocessed)

    output:
    tuple val(study_uid), path("results/"), emit: results

    script:
    """
    python -m maps.ct_head_hemorrhage.inference \
        --input-dir ${preprocessed} \
        --model-path ${params.model_dir}/hemorrhage-triage-v2.1.pt \
        --output-dir results/
    """
}

process POSTPROCESS {
    tag "$study_uid"
    label 'process_low'
    container 'imaging-agent/ct-head-hemorrhage:latest'

    input:
    tuple val(study_uid), path(results)

    output:
    tuple val(study_uid), path("output/"), emit: output

    script:
    """
    python -m maps.ct_head_hemorrhage.postprocess \
        --input-dir ${results} \
        --output-dir output/ \
        --postgres-url ${params.postgres_url} \
        --orthanc-url ${params.orthanc_url}
    """
}

workflow CT_HEAD_HEMORRHAGE {
    take:
    study_uid

    main:
    PREPROCESS(study_uid)
    INFERENCE(PREPROCESS.out.preprocessed)
    POSTPROCESS(INFERENCE.out.results)

    emit:
    output = POSTPROCESS.out.output
}
```

### nextflow.config

```groovy
params {
    // Defaults
    workflow         = null
    study_uid        = null
    input_dir        = null
    outdir           = "${launchDir}/results"

    // Service URLs
    orthanc_url      = 'http://localhost:8042'
    postgres_url     = 'postgresql://imaging:imaging_secret@localhost:5432/imaging_agent'
    nim_url          = 'http://localhost:8520/v1'

    // Model paths
    model_dir        = "${projectDir}/models"
}

profiles {
    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    dgx_spark {
        docker.enabled = true
        docker.runOptions = '--gpus all'

        process {
            withLabel: 'gpu' {
                containerOptions = '--gpus all'
                accelerator = 1
                memory = '64 GB'
            }
            withLabel: 'process_low' {
                cpus = 2
                memory = '4 GB'
            }
        }

        params {
            max_gpus = 1
            max_memory = '128.GB'
            max_cpus = 12
        }
    }
}

process {
    cpus = 4
    memory = '8 GB'
    time = '30m'

    errorStrategy = { task.exitStatus in [143, 137, 104, 134, 139] ? 'retry' : 'finish' }
    maxRetries = 2
}

report {
    enabled = true
    file = "${params.outdir}/pipeline_report.html"
}

timeline {
    enabled = true
    file = "${params.outdir}/timeline.html"
}

trace {
    enabled = true
    file = "${params.outdir}/trace.tsv"
}
```

---

## 19. Streamlit Portal

### services/portal/app.py

```python
"""Imaging Intelligence Agent — Streamlit Portal."""

import os

import requests
import streamlit as st

# ── Page Config ─────────────────────────────────────────────
st.set_page_config(
    page_title="Imaging Intelligence Agent",
    page_icon="🏥",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── NVIDIA Theme CSS ────────────────────────────────────────
st.markdown("""
<style>
    :root {
        --nvidia-green: #76b900;
        --dark-bg: #1a1a2e;
    }
    .main-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        padding: 2rem;
        border-radius: 16px;
        margin-bottom: 2rem;
        border-left: 5px solid #76b900;
    }
    .main-header h1 { color: white; margin: 0; font-size: 2rem; }
    .main-header p { color: #a0aec0; margin: 0.5rem 0 0 0; }
    .metric-card {
        background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%);
        color: white; padding: 1.5rem; border-radius: 12px; text-align: center;
    }
    .worklist-urgent { border-left: 4px solid #dc2626; }
    .worklist-moderate { border-left: 4px solid #f5a623; }
    .worklist-routine { border-left: 4px solid #059669; }
</style>
""", unsafe_allow_html=True)

# ── Service URLs ────────────────────────────────────────────
AGENT_URL = os.environ.get("AGENT_URL", "http://localhost:8524")
ORTHANC_URL = os.environ.get("ORTHANC_URL", "http://localhost:8042")
POSTGRES_URL = os.environ.get("POSTGRES_URL", "")


def get_service_host():
    """Detect host IP for external links."""
    host = os.environ.get("SERVICE_HOST")
    if host:
        return host
    try:
        import socket
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        host = s.getsockname()[0]
        s.close()
        return host
    except Exception:
        return "localhost"


SERVICE_HOST = get_service_host()


# ── Header ──────────────────────────────────────────────────
st.markdown("""
<div class="main-header">
    <h1>Imaging Intelligence Agent</h1>
    <p>CT / MRI / X-Ray — Automated Detection, Segmentation, Triage</p>
    <p style="color: #76b900; font-size: 0.9rem;">
        NVIDIA DGX Spark | MONAI Deploy | Open Architecture
    </p>
</div>
""", unsafe_allow_html=True)


# ── Sidebar ─────────────────────────────────────────────────
with st.sidebar:
    st.markdown("### Navigation")
    page = st.radio(
        "Select Page",
        ["Worklist", "Studies", "Agent Activity", "Monitoring"],
        label_visibility="collapsed",
    )

    st.markdown("### Service Status")
    services = {
        "Orthanc": f"{ORTHANC_URL}/system",
        "Agent": f"{AGENT_URL}/health",
    }
    for name, url in services.items():
        try:
            r = requests.get(url, timeout=2)
            st.markdown(f"🟢 {name}")
        except Exception:
            st.markdown(f"🔴 {name}")


# ── Pages ───────────────────────────────────────────────────
if page == "Worklist":
    st.markdown("## Active Worklist")
    # Query active_worklist view from PostgreSQL
    # Display as prioritized table with urgency color coding

elif page == "Studies":
    st.markdown("## Processed Studies")
    # Query study_summary view
    # Display with finding counts, status badges

elif page == "Agent Activity":
    st.markdown("## Agent Activity Log")
    # Query provenance table for recent runs
    # Show pipeline durations, model versions, status

elif page == "Monitoring":
    st.markdown("## GPU & Pipeline Monitoring")
    grafana_url = f"http://{SERVICE_HOST}:3000"
    st.markdown(f"[Open Grafana Dashboard]({grafana_url})")
    # Embed GPU metrics from DCGM exporter
```

---

## 20. Monitoring Stack

### config/prometheus.yml

```yaml
global:
  scrape_interval: 15s
  evaluation_interval: 15s
  external_labels:
    monitor: 'imaging-agent-monitor'

scrape_configs:
  - job_name: 'prometheus'
    static_configs:
      - targets: ['localhost:9090']

  - job_name: 'dcgm-exporter'
    static_configs:
      - targets: ['dcgm-exporter:9400']
        labels:
          instance: 'dgx-spark'
    scrape_interval: 5s

  - job_name: 'orthanc'
    static_configs:
      - targets: ['orthanc:8042']
    metrics_path: '/metrics'
    scrape_interval: 30s

  - job_name: 'embedding-service'
    static_configs:
      - targets: ['embedding-service:8000']
    metrics_path: '/metrics'

  - job_name: 'agent'
    static_configs:
      - targets: ['agent:8000']
    metrics_path: '/metrics'
```

### Key DCGM Metrics

| Metric | Description |
|---|---|
| `DCGM_FI_DEV_GPU_UTIL` | GPU utilization % |
| `DCGM_FI_DEV_GPU_TEMP` | GPU temperature (C) |
| `DCGM_FI_DEV_POWER_USAGE` | Power draw (W) |
| `DCGM_FI_DEV_FB_USED` | GPU memory used (MB) |
| `DCGM_FI_DEV_FB_FREE` | GPU memory free (MB) |

### Alert Rules (add to prometheus.yml)

```yaml
rule_files:
  - '/etc/prometheus/alerts.yml'

# alerts.yml
groups:
  - name: imaging-agent
    rules:
      - alert: HighGPUMemory
        expr: DCGM_FI_DEV_FB_USED / (DCGM_FI_DEV_FB_USED + DCGM_FI_DEV_FB_FREE) > 0.9
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "GPU memory usage above 90%"

      - alert: InferenceFailureRate
        expr: rate(inference_failures_total[5m]) > 0.1
        for: 2m
        labels:
          severity: critical
        annotations:
          summary: "High inference failure rate"
```

---

## 21. Testing Strategy

### Synthetic DICOM Generation

```python
"""Generate synthetic DICOM test data."""

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
from datetime import date


def create_synthetic_ct_head(output_path: str, num_slices: int = 64) -> str:
    """Generate a synthetic CT head DICOM series."""
    study_uid = generate_uid()
    series_uid = generate_uid()

    for i in range(num_slices):
        ds = Dataset()
        ds.PatientID = "TEST-001"
        ds.PatientName = "Test^Patient"
        ds.StudyDate = date.today().strftime("%Y%m%d")
        ds.Modality = "CT"
        ds.StudyDescription = "CT HEAD W/O CONTRAST"
        ds.BodyPartExamined = "HEAD"
        ds.StudyInstanceUID = study_uid
        ds.SeriesInstanceUID = series_uid
        ds.SOPInstanceUID = generate_uid()
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"  # CT Image Storage
        ds.Rows = 512
        ds.Columns = 512
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.RescaleIntercept = -1024
        ds.RescaleSlope = 1
        ds.PixelSpacing = [0.5, 0.5]
        ds.SliceThickness = 2.5
        ds.ImagePositionPatient = [0, 0, i * 2.5]
        ds.InstanceNumber = i + 1

        # Generate synthetic pixel data (noise + optional hemorrhage blob)
        pixels = np.random.randint(-100, 100, (512, 512), dtype=np.int16)
        ds.PixelData = pixels.tobytes()

        file_ds = FileDataset(
            f"{output_path}/slice_{i:04d}.dcm",
            ds,
            preamble=b"\x00" * 128,
        )
        file_ds.is_little_endian = True
        file_ds.is_implicit_VR = False
        file_ds.save_as(f"{output_path}/slice_{i:04d}.dcm")

    return study_uid
```

### pytest Fixtures

```python
"""Shared test fixtures — conftest.py."""

import pytest
import psycopg2
from testcontainers.postgres import PostgresContainer

from tests.synthetic_dicom import create_synthetic_ct_head


@pytest.fixture(scope="session")
def postgres_container():
    """Spin up PostgreSQL + pgvector for integration tests."""
    with PostgresContainer("pgvector/pgvector:pg16") as pg:
        conn = psycopg2.connect(pg.get_connection_url())
        with open("db/init.sql") as f:
            conn.cursor().execute(f.read())
        conn.commit()
        yield pg


@pytest.fixture
def db_conn(postgres_container):
    """Database connection for a single test."""
    conn = psycopg2.connect(postgres_container.get_connection_url())
    yield conn
    conn.rollback()
    conn.close()


@pytest.fixture
def synthetic_ct_head(tmp_path):
    """Generate synthetic CT head DICOM in a temp directory."""
    study_uid = create_synthetic_ct_head(str(tmp_path / "ct_head"), num_slices=32)
    return {"path": str(tmp_path / "ct_head"), "study_uid": study_uid}
```

### Example Test

```python
"""Test CT Head hemorrhage workflow."""

import numpy as np
from maps.ct_head_hemorrhage.operators import classify_urgency, measure_midline_shift


def test_urgency_critical():
    assert classify_urgency(volume_ml=35.0, shift_mm=6.0, thickness_mm=12.0) == ("critical", "P1")


def test_urgency_urgent():
    assert classify_urgency(volume_ml=12.0, shift_mm=3.0, thickness_mm=8.0) == ("urgent", "P2")


def test_urgency_routine():
    assert classify_urgency(volume_ml=2.0, shift_mm=1.0, thickness_mm=3.0) == ("routine", "P4")


def test_midline_shift_no_hemorrhage():
    seg = np.zeros((64, 64, 64), dtype=np.uint8)
    result = measure_midline_shift(seg, voxel_spacing=(1.0, 1.0, 1.0))
    assert result["shift_mm"] == 0.0


def test_midline_shift_left():
    seg = np.zeros((64, 64, 64), dtype=np.uint8)
    seg[10:20, 30:40, 30:40] = 1  # Left-sided hemorrhage
    result = measure_midline_shift(seg, voxel_spacing=(1.0, 1.0, 1.0))
    assert result["shift_mm"] > 0
    assert result["shift_direction"] == "leftward"
```

---

## 22. ARM64 Compatibility Guide

### Base Images (ARM64-compatible)

| Image | ARM64 Support | Notes |
|---|---|---|
| `nvcr.io/nvidia/pytorch:24.01-py3` | Yes | Check NGC for `-aarch64` tags |
| `pgvector/pgvector:pg16` | Yes | Multi-arch |
| `orthancteam/orthanc:24.1.2` | Yes | Multi-arch |
| `prom/prometheus:v2.48.0` | Yes | Multi-arch |
| `grafana/grafana:10.2.2` | Yes | Multi-arch |
| `python:3.11-slim` | Yes | Multi-arch |

### Python Packages with ARM64 Wheels

All of these have pre-built ARM64 wheels on PyPI:
- `torch` (via NGC container or `pip install torch`)
- `monai` (pure Python + optional C extensions)
- `numpy`, `scipy`, `scikit-image`
- `pydicom`, `highdicom`
- `psycopg2-binary` (use `psycopg2-binary` not `psycopg2`)
- `langchain`, `langgraph` (pure Python)
- `fastapi`, `uvicorn`
- `streamlit`
- `transformers`
- `pgvector`

### NIM Container Variants

For DGX Spark, append `-dgx-spark` to standard NIM image tags:

```bash
# Standard (x86_64)
nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest

# DGX Spark (ARM64)
nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark
```

### Docker Buildx for Multi-arch

```bash
# Build ARM64 images on x86_64 dev machine
docker buildx create --use --name multiarch
docker buildx build --platform linux/arm64 -t imaging-agent/ct-head:latest --push .
```

---

## 23. Configuration Reference

### .env.example

```bash
# ── PostgreSQL ──────────────────────────────────────────────
POSTGRES_USER=imaging
POSTGRES_PASSWORD=imaging_secret
POSTGRES_DB=imaging_agent

# ── NVIDIA ──────────────────────────────────────────────────
NGC_API_KEY=your_ngc_api_key_here

# ── Monitoring ──────────────────────────────────────────────
GRAFANA_USER=admin
GRAFANA_PASSWORD=changeme

# ── FHIR (optional) ────────────────────────────────────────
FHIR_SERVER_URL=http://localhost:8080/fhir

# ── Service Host (auto-detected if not set) ─────────────────
# SERVICE_HOST=192.168.1.100
```

---

## 24. Implementation Sequence

Build in this order — each step depends on the ones above it:

| Step | Component | Dependencies | Verification |
|---|---|---|---|
| 1 | PostgreSQL + pgvector | None | `psql` connect, `\dx` shows vector extension |
| 2 | Orthanc DICOM server | None | `curl http://localhost:8042/system` returns JSON |
| 3 | DICOM listener service | Orthanc, PostgreSQL | Upload DICOM → webhook fires → study row in DB |
| 4 | CXR Rapid Findings MAP | Orthanc, PostgreSQL | Upload CXR → findings in DB, < 30s |
| 5 | CT Head Hemorrhage MAP | Orthanc, PostgreSQL | Upload CT head → findings + measurements in DB |
| 6 | Embedding service | PostgreSQL (pgvector) | POST /embed → 384-dim vector stored |
| 7 | LangGraph agent | PostgreSQL, NIM LLM | Query findings, generate summary |
| 8 | FHIR publisher | PostgreSQL | POST /publish → FHIR Bundle created |
| 9 | Streamlit portal | PostgreSQL, Orthanc, Agent | Worklist renders, study list populates |
| 10 | CT Chest Lung Nodule MAP | Orthanc, PostgreSQL, priors | Upload CT → nodule tracking with VDT |
| 11 | MRI Brain MS Lesion MAP | Orthanc, PostgreSQL, priors | Upload MRI → lesion tracking with registration |
| 12 | NIM LLM + RAG pipeline | NIM container, pgvector | RAG query returns grounded response |
| 13 | Nextflow orchestration | All MAPs, all services | `nextflow run main.nf -profile dgx_spark --workflow ct_head` |
| 14 | Monitoring stack | Prometheus, Grafana, DCGM | Dashboard shows GPU metrics and pipeline traces |
| 15 | End-to-end integration | Everything | Upload 4 studies → all findings + SR + FHIR + worklist |

### Quick Start Commands

```bash
# 1. Clone and configure
git clone <repo-url> hls-imaging-agent
cd hls-imaging-agent
cp .env.example .env
# Edit .env with your NGC_API_KEY

# 2. Start infrastructure
docker compose up -d orthanc postgres

# 3. Download models
bash models/download_models.sh

# 4. Start all services
docker compose up -d

# 5. Upload test DICOM
bash scripts/seed_test_data.sh

# 6. Run a workflow
nextflow run main.nf -profile dgx_spark --workflow cxr --study_uid <uid>

# 7. Open portal
# http://localhost:8525
```

---

## Appendix A — Technology Stack Summary

| Component | Role | License |
|---|---|---|
| MONAI | Medical imaging AI framework | Apache 2.0 |
| MONAI Deploy | Containerized inference packaging (MAPs) | Apache 2.0 |
| NVIDIA FLARE | Federated learning across institutions | Apache 2.0 |
| Nextflow | Pipeline DAG orchestration | Apache 2.0 |
| LangChain / LangGraph | Agent orchestration + MCP tools | MIT |
| PostgreSQL + pgvector | Structured + vector query | PostgreSQL License |
| Orthanc | DICOM server | GPLv3 |
| dcm4chee | DICOM archive (alternative) | MPL 1.1 / GPL 2.0 / LGPL 2.1 |
| pydicom | DICOM parsing and manipulation | MIT |
| highdicom | DICOM SR / SEG construction | MIT |
| HAPI FHIR | FHIR server for clinical integration | Apache 2.0 |
| fhir.resources | Pydantic FHIR models | BSD |
| Grafana + Prometheus | Monitoring and observability | AGPL / Apache 2.0 |
| NVIDIA NIM | Inference microservices | NVAIE ($4,500/GPU/yr) |
| NVIDIA Parabricks | GPU-accelerated genomics | NVAIE ($4,500/GPU/yr) |
| NVIDIA BioNeMo | Drug discovery + molecular modeling | NVAIE ($4,500/GPU/yr) |

---

*HCLS AI Factory — Open Source (Apache 2.0) — NVIDIA DGX Spark*
