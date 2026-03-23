# Single-Cell Intelligence Agent -- Documentation Index

> **Part of the [Precision Intelligence Network](../engines/precision-intelligence.md)** — one of 11 specialized agents sharing a common molecular foundation within the HCLS AI Factory.

**Version:** 1.0.0
**Date:** 2026-03-22
**Agent:** Single-Cell Intelligence Agent
**Ports:** API 8540 | UI 8130

---

## Documentation Set

| # | Document | Description | Audience |
|---|---------|-------------|----------|
| 1 | [Production Readiness Report](PRODUCTION_READINESS_REPORT.md) | 25-section PRR covering architecture, collections, workflows, decision engines, security, testing, and go/no-go recommendation | Engineering, QA |
| 2 | [Project Bible](PROJECT_BIBLE.md) | Mission, problem statement, solution architecture, feature inventory, technical stack, data flow, quality gates, and roadmap | All stakeholders |
| 3 | [Architecture Guide](ARCHITECTURE_GUIDE.md) | Layered architecture, GPU acceleration pipeline, RAPIDS integration, TME classification pipeline, spatial deconvolution, RAG search architecture | Engineering |
| 4 | [White Paper](WHITE_PAPER.md) | Resolution gap in precision medicine, tumor heterogeneity, GPU necessity, clinical applications, foundation models | Executive, clinical |
| 5 | [Deployment Guide](DEPLOYMENT_GUIDE.md) | Standalone Docker Compose, integrated DGX Spark, manual development setup, configuration reference, monitoring, scaling | DevOps, engineering |
| 6 | [Demo Guide](DEMO_GUIDE.md) | 5 demos: cell type annotation, TME profiling, drug response, CAR-T validation, spatial analysis | Sales, clinical, demos |
| 7 | [Learning Guide: Foundations](LEARNING_GUIDE_FOUNDATIONS.md) | Single-cell primer: technologies, data formats (AnnData), analysis pipeline (QC to clustering to DE) | New team members |
| 8 | [Learning Guide: Advanced](LEARNING_GUIDE_ADVANCED.md) | TME classification deep dive, subclonal architecture, spatial transcriptomics, trajectory inference, foundation models (scGPT/Geneformer), GPU benchmarks | Experienced analysts |
| 9 | [Research Paper](SINGLE_CELL_INTELLIGENCE_AGENT_RESEARCH_PAPER.md) | Academic research paper on the agent's design and capabilities | Publications |

---

## Quick Reference

### Agent Statistics

| Metric | Value |
|--------|-------|
| Milvus collections | 12 |
| Analysis workflows | 10 |
| Decision support engines | 4 |
| Cell types in knowledge base | 44 |
| Drugs modeled | 30 |
| Marker genes | 75 |
| Immune signatures | 10 |
| Ligand-receptor pairs | 25 |
| Cancer TME atlas profiles | 12 |
| API endpoints | 25 |
| Test cases | ~185 |
| Source code lines | 14,560 |

### Port Map

| Port | Service |
|------|---------|
| 8540 | FastAPI REST API |
| 8130 | Streamlit UI |
| 19530 | Milvus (shared) |
| 69530 | Milvus (standalone) |

### Key Files

| Path | Purpose |
|------|---------|
| `config/settings.py` | All configuration with `SC_` env prefix |
| `src/agent.py` | Core reasoning engine |
| `src/collections.py` | 12 Milvus collection schemas |
| `src/clinical_workflows.py` | 10 analysis workflows |
| `src/decision_support.py` | 4 clinical engines |
| `src/knowledge.py` | Domain knowledge base |
| `api/main.py` | FastAPI application |
| `app/sc_ui.py` | Streamlit UI |
| `docker-compose.yml` | Standalone deployment |

---

*HCLS AI Factory -- Single-Cell Intelligence Agent Documentation Index v1.0.0*

---

!!! warning "Clinical Decision Support Disclaimer"
    This agent is a clinical decision support research tool. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
