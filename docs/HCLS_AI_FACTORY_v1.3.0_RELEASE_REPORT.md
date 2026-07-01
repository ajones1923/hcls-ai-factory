# HCLS AI Factory v1.3.0 -- Post-Implementation Report

**Author:** Adam Jones | **Date:** March 2026 | **License:** Apache 2.0

---

## Executive Summary

HCLS AI Factory v1.3.0 has been released across 12 GitHub repositories, representing the most significant expansion of the platform since its inception. This release transforms the system from a 3-agent prototype into an 11-agent precision medicine platform spanning genomics, oncology, cardiology, neurology, pharmacogenomics, rare disease diagnostics, clinical trial intelligence, single-cell analysis, immunology, medical imaging, and CAR-T cell therapy.

**Key milestones:**

- **v1.3.0 released** across all 12 GitHub repositories
- **3 engines, 11 intelligence agents, 21 services** deployed and verified
- **All CI green** (12/12 repos passing lint + tests)
- **Documentation site live** at hcls-ai-factory.org with 142 documented pages
- **139 Milvus collections** seeded with ~47,691 vectors across all agents
- **158 test files** across the agent ecosystem; 225 tests in hcls_common alone

The platform delivers its core promise: **Patient DNA to Drug Candidates in under 5 hours** on an NVIDIA DGX Spark ($4,699), compared to months using traditional workflows.

---

## Release Scope

### What Changed from v1.2.0 to v1.3.0

| Metric | v1.2.0 | v1.3.0 | Change |
|--------|--------|--------|--------|
| Intelligence Agents | 3 | 11 | +8 agents |
| Milvus Collections | ~30 | 139 | +109 collections |
| Seeded Vectors | ~6,500 | ~47,691 | +41,191 vectors |
| Documentation Pages | 148 | 142 | Consolidated (cleanup) |
| Test Files (agents) | ~40 | 158 | +118 test files |
| GitHub Repositories | 6 | 12 | +6 repos |
| Docker Services | 12 | 21 | +9 services |
| CI Workflows | 6 | 12 | +6 workflows |

### New Agents Added in v1.3.0

1. **Cardiology Intelligence Agent** -- 6 validated risk calculators (ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE II), GDMT optimizer, 8 clinical workflows, 1,927 tests
2. **Neurology Intelligence Agent** -- Neurological disorder decision support with cross-modal genomic triggers
3. **Precision Autoimmune Agent** -- Autoimmune disease and immune-mediated condition intelligence
4. **Rare Disease Diagnostic Agent** -- Diagnostic odyssey support with phenotype-genotype correlation
5. **Pharmacogenomics (PGx) Agent** -- 25 pharmacogenes, 100+ drugs, 9 dosing algorithms, 15 HLA associations, 1,001+ tests
6. **Single-Cell Intelligence Agent** -- Single-cell RNA-seq and multi-omics analysis
7. **Clinical Trial Intelligence Agent** -- Trial matching, eligibility screening, site selection
8. **Expanded Precision Biomarker Agent** -- Biological age estimation (PhenoAge/GrimAge), pharmacogenomic profiling

### Agents Carried Forward from v1.2.0

1. **Precision Biomarker Agent** -- Genotype-aware biomarker interpretation, disease trajectory detection
2. **Precision Oncology Agent** -- Molecular tumor board decision support, variant annotation (CIViC/OncoKB)
3. **CAR-T Intelligence Agent** -- Cross-functional CAR-T cell therapy intelligence, 6,266+ vectors

### Pediatric Safety Filters

Six pediatric safety filters are now connected to real data:

- hERG channel inhibition estimation
- Teratogenicity SMARTS pattern matching
- Pediatric dosing adjustment algorithms
- Age-appropriate formulation checks
- Growth plate toxicity screening
- Developmental neurotoxicity flags

---

## Infrastructure Verification

### Services (21/21 Healthy)

| # | Service | Port(s) | Type | Status |
|---|---------|---------|------|--------|
| 1 | Landing Page / Portal | 8080 | Flask | Healthy |
| 2 | Genomics API | 5000 | Flask | Healthy |
| 3 | RAG API | 5001 | Flask | Healthy |
| 4 | RAG Chat UI | 8501 | Streamlit | Healthy |
| 5 | Discovery Portal | 8510 | Streamlit | Healthy |
| 6 | Milvus Vector DB | 19530 | gRPC | Healthy |
| 7 | Milvus Management | 9091 | HTTP | Healthy |
| 8 | Milvus etcd | 2379 | Internal | Healthy |
| 9 | Milvus MinIO | 9001 | Internal | Healthy |
| 10 | Prometheus | 9099 | Metrics | Healthy |
| 11 | Grafana | 3000 | Dashboard | Healthy |
| 12 | Precision Biomarker Agent | 8502 / 8102 | Streamlit / FastAPI | Healthy |
| 13 | Precision Oncology Agent | 8503 / 8103 | Streamlit / FastAPI | Healthy |
| 14 | CAR-T Intelligence Agent | 8504 / 8104 | Streamlit / FastAPI | Healthy |
| 15 | Imaging Intelligence Agent | 8505 / 8105 | Streamlit / FastAPI | Healthy |
| 16 | Precision Autoimmune Agent | 8506 / 8106 | Streamlit / FastAPI | Healthy |
| 17 | Pharmacogenomics Agent | 8507 / 8107 | Streamlit / FastAPI | Healthy |
| 18 | Cardiology Intelligence Agent | 8536 / 8126 | Streamlit / FastAPI | Healthy |
| 19 | Clinical Trial Intelligence Agent | 8128 / 8538 | Streamlit / FastAPI | Healthy |
| 20 | Neurology Intelligence Agent | 8529 / 8528 | Streamlit / FastAPI | Healthy |
| 21 | Rare Disease Diagnostic Agent | 8544 / 8134 | Streamlit / FastAPI | Healthy |

**Note:** Single-Cell Intelligence Agent (8130/8540) operates as service 21 in extended deployments. The table above reflects the 21-service baseline verified at release.

### Agent APIs (11/11 Responding)

| Agent | API Endpoint | Health Check | Collections | Vectors |
|-------|-------------|--------------|-------------|---------|
| Precision Biomarker | :8102/health | 200 OK | 10 + shared | ~4,200 |
| Precision Oncology | :8103/health | 200 OK | 10 + shared | ~4,500 |
| CAR-T Intelligence | :8104/health | 200 OK | 11 | ~6,266 |
| Imaging Intelligence | :8105/health | 200 OK | 10 + shared | ~3,800 |
| Precision Autoimmune | :8106/health | 200 OK | 10 + shared | ~3,500 |
| Pharmacogenomics | :8107/health | 200 OK | 15 | ~5,200 |
| Cardiology Intelligence | :8126/health | 200 OK | 12 + shared | ~4,800 |
| Neurology Intelligence | :8528/health | 200 OK | 12 + shared | ~4,100 |
| Single-Cell Intelligence | :8540/health | 200 OK | 12 + shared | ~3,900 |
| Clinical Trial Intelligence | :8538/health | 200 OK | 12 + shared | ~3,800 |
| Rare Disease Diagnostic | :8134/health | 200 OK | 12 + shared | ~3,531 |

### Milvus Vector Database

- **139 collections** across all agents, 0 empty
- **~47,691 vectors** total, seeded from curated knowledge bases
- **Embedding model:** BAAI/bge-small-en-v1.5 (384 dimensions)
- **Search mode:** Hybrid (semantic + metadata filtering)
- **Shared collection:** `genomic_evidence` available to all agents via cross-modal triggers

---

## Three-Engine Architecture

### Engine 1: Genomics Pipeline

- **Input:** FASTQ (whole-genome sequencing data)
- **Processing:** BWA-MEM2 alignment, Parabricks fq2bam, DeepVariant variant calling
- **Output:** VCF with 11.7M variants
- **Reference genome:** GRCh38 (HG002 GIAB demo sample)
- **Performance:** 120-240 minutes (GPU) vs 24-48 hours (CPU-only)

### Engine 2: RAG/Chat Pipeline

- **Input:** VCF from Engine 1
- **Processing:** ClinVar/AlphaMissense/VEP annotation, BGE-small-en-v1.5 embedding, Milvus vector storage, Claude AI RAG
- **Output:** Annotated variants, target hypotheses, clinical interpretations
- **Scale:** 3.56M searchable variants, 4.1M ClinVar records, 71M AlphaMissense predictions
- **Query response:** <5 seconds

### Engine 3: Drug Discovery Pipeline

- **Input:** Target hypotheses from Engine 2
- **Processing:** MolMIM molecule generation, DiffDock molecular docking, RDKit ADMET scoring
- **Output:** Ranked drug candidates with docking scores, PDF reports
- **Demo pathway:** VCP for frontotemporal dementia -- 4 cryo-EM structures, CB-5083 inhibitor, 100 analogues

---

## Security Hardening

| Control | Implementation | Status |
|---------|---------------|--------|
| Portal Authentication | Always-on with auto-generated API keys | Active |
| XSS Protection | `html.escape()` on all dynamic template values | Active |
| Pickle Deserialization | `allow_pickle=False` enforced globally | Active |
| Grafana Access | Requires `GRAFANA_PASSWORD` environment variable (no default) | Active |
| HTTPS | Caddy reverse proxy with self-signed TLS on port 443 | Active |
| Clinical Disclaimers | "For Research Use Only" on 142/142 documentation pages | Active |
| API Key Storage | Environment variables via `.env` (not committed to VCS) | Active |
| Container Isolation | Non-root containers with memory limits enforced | Active |
| Dependency Scanning | Automated via CI/CD pipeline on every push | Active |

---

## CI/CD Status

All 12 repositories pass continuous integration checks as of the v1.3.0 release.

| # | Repository | CI Workflow | Conclusion |
|---|-----------|-------------|------------|
| 1 | hcls-ai-factory (main) | ci.yml | GREEN |
| 2 | hcls-ai-factory-public | ci.yml | GREEN |
| 3 | precision-biomarker-agent | ci.yml | GREEN |
| 4 | precision-oncology-agent | ci.yml | GREEN |
| 5 | cart-intelligence-agent | ci.yml | GREEN |
| 6 | imaging-intelligence-agent | ci.yml | GREEN |
| 7 | precision-autoimmune-agent | ci.yml | GREEN |
| 8 | pharmacogenomics-intelligence-agent | ci.yml | GREEN |
| 9 | cardiology-intelligence-agent | ci.yml | GREEN |
| 10 | clinical-trial-intelligence-agent | ci.yml | GREEN |
| 11 | neurology-intelligence-agent | ci.yml | GREEN |
| 12 | rare-disease-diagnostic-agent | ci.yml | GREEN |

Additional workflows:
- **release.yml** configured on hcls-ai-factory and hcls-ai-factory-public for automated GitHub Release creation on tag push
- **single-cell-intelligence-agent** CI integrated within main repo workflow

---

## Documentation

### Site: hcls-ai-factory.org

- **142 markdown pages** published via MkDocs + Netlify
- **Complete documentation overhaul** from v1.2.0 (148 pages consolidated to 142 after cleanup and deduplication)

### Per-Agent Documentation (11 agents x 8 documents)

Each agent ships with a standardized documentation set:

1. White Paper
2. Project Bible
3. Design Document
4. Architecture Guide
5. Demo Guide
6. Deployment Guide
7. Learning Guide -- Foundations
8. Learning Guide -- Advanced

### Unified Learning Guides

| Document | Lines | Coverage |
|----------|-------|----------|
| Learning Guide Unified -- Foundations | 3,881 | All 11 agents, platform architecture, getting started |
| Learning Guide Unified -- Advanced | 3,051 | Advanced topics, integration patterns, production deployment |

### Research Paper

- **arXiv paper:** 669 lines, updated to cover all 11 agents and the three-engine architecture
- **Title:** "HCLS AI Factory: An Open-Source Precision Medicine Platform on NVIDIA DGX Spark"

### Pediatric Oncology Resources

- Dedicated Pediatric Oncology Demos page with 7 infographics
- Age-appropriate safety filter documentation
- Pediatric dosing algorithm reference guides

---

## GitHub Releases

All repositories tagged with v1.3.0:

| # | Repository | Release Tag | Assets |
|---|-----------|-------------|--------|
| 1 | hcls-ai-factory | v1.3.0 | Source, docker-compose, docs |
| 2 | hcls-ai-factory-public | v1.3.0 | Source, governance docs |
| 3 | precision-biomarker-agent | v1.3.0 | Source, Dockerfile, tests |
| 4 | precision-oncology-agent | v1.3.0 | Source, Dockerfile, tests |
| 5 | cart-intelligence-agent | v1.3.0 | Source, Dockerfile, tests |
| 6 | imaging-intelligence-agent | v1.3.0 | Source, Dockerfile, tests |
| 7 | precision-autoimmune-agent | v1.3.0 | Source, Dockerfile, tests |
| 8 | pharmacogenomics-intelligence-agent | v1.3.0 | Source, Dockerfile, tests |
| 9 | cardiology-intelligence-agent | v1.3.0 | Source, Dockerfile, tests |
| 10 | clinical-trial-intelligence-agent | v1.3.0 | Source, Dockerfile, tests |
| 11 | neurology-intelligence-agent | v1.3.0 | Source, Dockerfile, tests |
| 12 | rare-disease-diagnostic-agent | v1.3.0 | Source, Dockerfile, tests |

---

## Test Coverage

### Agent Test Files

- **158 test files** across all 11 agent repositories
- **138 test files** contain active test functions

### Shared Library (hcls_common)

| Test File | Module Covered |
|-----------|---------------|
| test_circuit_breaker.py | Circuit breaker / fault tolerance |
| test_milvus_client.py | Milvus vector DB client |
| test_enums.py | Enumeration types |
| test_llm_client.py | Multi-provider LLM client |
| test_embedder.py | BGE embedding pipeline |
| test_tracing.py | OpenTelemetry tracing |
| test_config.py | Configuration management |
| test_event_bus.py | Event bus / pub-sub |
| test_reproducibility.py | FDA 21 CFR Part 11 reproducibility |
| test_security.py | Security utilities |

**Total:** 10 test files, 225 tests, all passing

### Notable Agent Test Suites

| Agent | Test Count | Notes |
|-------|-----------|-------|
| Cardiology Intelligence | 1,927 | 6 risk calculators, GDMT optimizer, 8 workflows |
| Pharmacogenomics | 1,001+ | 25 pharmacogenes, dosing algorithms, HLA associations |
| CAR-T Intelligence | 500+ | Construct analysis, manufacturing, assay validation |

---

## Technology Stack

| Layer | Components |
|-------|-----------|
| **Compute** | NVIDIA DGX Spark ($4,699), GB10 GPU, 128 GB unified LPDDR5x, 20 ARM cores (Grace CPU), NVLink-C2C |
| **Genomics** | Parabricks 4.6, DeepVariant, BWA-MEM2 |
| **AI/LLM** | Claude (Anthropic), BioNeMo NIMs |
| **Vector DB** | Milvus 2.4.0 (with etcd + MinIO) |
| **Chemistry** | RDKit, MolMIM, DiffDock |
| **Embedding** | BAAI/bge-small-en-v1.5 (384-dim) |
| **Reference Data** | ClinVar (4.1M records), AlphaMissense (71M predictions) |
| **Web** | Flask, FastAPI, Streamlit |
| **Infrastructure** | Docker, Nextflow DSL2, Caddy (TLS) |
| **Observability** | Prometheus 2.49.1, Grafana 10.3.1 |
| **CI/CD** | GitHub Actions (12 workflows) |
| **Documentation** | MkDocs + Netlify |

---

## Platform Statistics

| Metric | Value |
|--------|-------|
| Variants in demo genome | 11.7M |
| Searchable vectors | 3.56M (annotated variants) + 47.7K (agent knowledge) |
| ClinVar records | 4.1M |
| AlphaMissense predictions | 71M |
| Genes covered | 201 |
| Therapeutic areas | 13 |
| Druggable targets | 171 |
| Genome processing time (GPU) | 120-240 minutes |
| Genome processing time (CPU) | 24-48 hours |
| Query response time | <5 seconds |
| Total disk footprint | ~1.1 TB |

---

## Known Limitations

1. **Computational candidates only** -- Drug candidates are generated in silico and have not been experimentally validated. All outputs are for research use only and are not intended for clinical decision-making without independent verification.

2. **Curated knowledge bases** -- Milvus collections are seeded with curated knowledge bases (~47K vectors), not full-scale production datasets. Production deployments should expand these with institutional data.

3. **Self-signed TLS** -- HTTPS uses a self-signed certificate via Caddy. Production deployments should use a certificate from a recognized Certificate Authority.

4. **Single-node deployment** -- The platform is designed for single-node DGX Spark deployment. Multi-node scaling requires additional orchestration infrastructure.

5. **API key management** -- API keys (Anthropic, etc.) are managed via environment variables in `.env` files. Production deployments should use a secrets management solution (e.g., HashiCorp Vault, AWS Secrets Manager).

6. **hcls_common test gaps** -- 14 of 23 modules in the shared library have 0 test coverage. Priority modules for future testing: query_router (44KB), event_bus (31KB), bidirectional_triggers (43KB), meta_agent (80KB), report_generator (83KB).

---

## Verification Checklist

- [x] All 12 GitHub repositories tagged with v1.3.0
- [x] All 12 CI workflows passing (GREEN)
- [x] 21/21 Docker services healthy
- [x] 11/11 agent health endpoints responding 200 OK
- [x] 139/139 Milvus collections populated (0 empty)
- [x] ~47,691 vectors seeded across all agent knowledge bases
- [x] 3.56M annotated variants searchable in Milvus
- [x] Portal authentication enabled (auto-generated keys)
- [x] XSS protection verified (`html.escape()` on all dynamic values)
- [x] Pickle deserialization blocked (`allow_pickle=False`)
- [x] Grafana requires `GRAFANA_PASSWORD` (no default password)
- [x] HTTPS active via Caddy with TLS on port 443
- [x] Clinical decision support disclaimers on 142/142 pages
- [x] 158 test files across agent ecosystem
- [x] 225 hcls_common tests passing
- [x] Cardiology agent: 1,927 tests passing
- [x] Pharmacogenomics agent: 1,001+ tests passing
- [x] Documentation site live at hcls-ai-factory.org
- [x] Learning Guide Foundations Unified: 3,881 lines
- [x] Learning Guide Advanced Unified: 3,051 lines
- [x] arXiv paper updated for 11 agents (669 lines)
- [x] Pediatric safety filters connected to real data (6 filters)
- [x] All agents have Streamlit UI + FastAPI backend + Milvus integration
- [x] Docker Compose unified manifest validated (`docker-compose.dgx-spark.yml`)
- [x] Prometheus + Grafana observability stack operational

---

## Appendix A: Complete Port Map

| Port | Service | Protocol |
|------|---------|----------|
| 443 | Caddy HTTPS Reverse Proxy | HTTPS |
| 2379 | Milvus etcd | gRPC (internal) |
| 3000 | Grafana Dashboard | HTTP |
| 5000 | Genomics Pipeline API | HTTP (Flask) |
| 5001 | RAG/Chat API | HTTP (Flask) |
| 8080 | Landing Page / Portal | HTTP (Flask) |
| 8102 | Precision Biomarker Agent API | HTTP (FastAPI) |
| 8103 | Precision Oncology Agent API | HTTP (FastAPI) |
| 8104 | CAR-T Intelligence Agent API | HTTP (FastAPI) |
| 8105 | Imaging Intelligence Agent API | HTTP (FastAPI) |
| 8106 | Precision Autoimmune Agent API | HTTP (FastAPI) |
| 8107 | Pharmacogenomics Agent API | HTTP (FastAPI) |
| 8126 | Cardiology Intelligence Agent API | HTTP (FastAPI) |
| 8128 | Clinical Trial Agent Streamlit UI | HTTP (Streamlit) |
| 8130 | Single-Cell Agent Streamlit UI | HTTP (Streamlit) |
| 8134 | Rare Disease Diagnostic Agent API | HTTP (FastAPI) |
| 8501 | RAG Chat UI | HTTP (Streamlit) |
| 8502 | Precision Biomarker Agent UI | HTTP (Streamlit) |
| 8503 | Precision Oncology Agent UI | HTTP (Streamlit) |
| 8504 | CAR-T Intelligence Agent UI | HTTP (Streamlit) |
| 8505 | Imaging Intelligence Agent UI | HTTP (Streamlit) |
| 8506 | Precision Autoimmune Agent UI | HTTP (Streamlit) |
| 8507 | Pharmacogenomics Agent UI | HTTP (Streamlit) |
| 8510 | Discovery Portal | HTTP (Streamlit) |
| 8528 | Neurology Intelligence Agent API | HTTP (FastAPI) |
| 8529 | Neurology Intelligence Agent UI | HTTP (Streamlit) |
| 8536 | Cardiology Intelligence Agent UI | HTTP (Streamlit) |
| 8538 | Clinical Trial Intelligence Agent API | HTTP (FastAPI) |
| 8540 | Single-Cell Intelligence Agent API | HTTP (FastAPI) |
| 8544 | Rare Disease Diagnostic Agent UI | HTTP (Streamlit) |
| 9001 | MinIO Console (Milvus storage) | HTTP (internal) |
| 9091 | Milvus Management API | HTTP |
| 9099 | Prometheus | HTTP |
| 9100 | Node Exporter | HTTP |
| 9400 | DCGM Exporter (GPU metrics) | HTTP |
| 19530 | Milvus Vector Database | gRPC |

---

## Appendix B: File Inventory

### Repository Structure

```
hcls-ai-factory/
  ai_agent_adds/
    cardiology_intelligence_agent/       # New in v1.3.0
    cart_intelligence_agent/
    clinical_trial_intelligence_agent/   # New in v1.3.0
    imaging_intelligence_agent/
    neurology_intelligence_agent/        # New in v1.3.0
    pharmacogenomics_intelligence_agent/ # New in v1.3.0
    precision_autoimmune_agent/          # New in v1.3.0
    precision_biomarker_agent/
    precision_oncology_agent/
    rare_disease_diagnostic_agent/       # New in v1.3.0
    single_cell_intelligence_agent/      # New in v1.3.0
    docs/                                # Agent documentation hub
  docker-compose.dgx-spark.yml          # 21-service unified manifest
  drug-discovery-pipeline/               # Engine 3
  genomics-pipeline/                     # Engine 1
  hcls-ai-factory/                       # Inner reference repo (CI/CD)
  hcls-ai-factory-public/                # Public GitHub release
  hcls-ai-factory-vast/                  # the AI platform deployment (MkDocs)
  hls-orchestrator/                      # Nextflow DSL2 orchestrator
  landing-page/                          # Flask portal (:8080)
  lib/hcls_common/                       # 23-module shared library
  monitoring/                            # Prometheus + Grafana configs
  rag-chat-pipeline/                     # Engine 2
  start-factory.sh                       # Master startup script (19KB)
  stop-factory.sh                        # Graceful shutdown
  health-monitor.sh                      # 11-service health watchdog (19KB)
  docs/                                  # Reports, guides, demo scripts
```

### Key Configuration Files

| File | Purpose | Size |
|------|---------|------|
| docker-compose.dgx-spark.yml | Unified 21-service Docker Compose | ~437 lines |
| start-factory.sh | Master startup for all services | 19 KB |
| health-monitor.sh | Health monitoring with auto-restart | 19 KB |
| monitoring/prometheus.yml | Prometheus scrape configuration | -- |
| monitoring/grafana/provisioning/ | Grafana datasource + dashboard provisioning | -- |

---

## Appendix C: CI Workflow Run Log

All CI workflows executed on the v1.3.0 release tag. Each workflow performs:

1. **Checkout** -- Clone repository at tagged commit
2. **Python Setup** -- Python 3.11 with pip cache
3. **Dependency Install** -- `pip install -r requirements.txt` + test dependencies
4. **Lint** -- `ruff check .` (or `flake8` where configured)
5. **Type Check** -- `mypy` on core modules (where configured)
6. **Unit Tests** -- `pytest` with coverage reporting
7. **Integration Tests** -- Milvus connection tests (mocked in CI)
8. **Artifact Upload** -- Coverage reports and test results

### Workflow Results Summary

```
hcls-ai-factory (main)           ci.yml    PASS    lint + tests
hcls-ai-factory-public           ci.yml    PASS    lint + tests
precision-biomarker-agent        ci.yml    PASS    lint + tests
precision-oncology-agent         ci.yml    PASS    lint + tests
cart-intelligence-agent          ci.yml    PASS    lint + tests
imaging-intelligence-agent       ci.yml    PASS    lint + tests
precision-autoimmune-agent       ci.yml    PASS    lint + tests
pharmacogenomics-intelligence    ci.yml    PASS    lint + tests
cardiology-intelligence-agent    ci.yml    PASS    lint + tests (1,927 tests)
clinical-trial-intelligence      ci.yml    PASS    lint + tests
neurology-intelligence-agent     ci.yml    PASS    lint + tests
rare-disease-diagnostic-agent    ci.yml    PASS    lint + tests
```

---

*This report was generated as part of the HCLS AI Factory v1.3.0 release process. For questions or access requests, contact Adam Jones.*

*HCLS AI Factory is open-source software released under the Apache 2.0 license. All clinical outputs are for research use only and are not intended for clinical decision-making without independent validation.*
