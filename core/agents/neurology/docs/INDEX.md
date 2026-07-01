# Neurology Intelligence Agent -- Documentation Index

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones

---

## Overview

The Neurology Intelligence Agent is a RAG-powered clinical decision support system covering 10 neurological disease domains, 14 Milvus vector collections, 10 validated clinical scale calculators, and 8 evidence-based clinical workflows. It is part of the HCLS AI Factory precision medicine platform.

**Ports:** FastAPI API (8528) | Streamlit UI (8529)

---

## Documentation Set

| # | Document | Description | Audience |
|---|---|---|---|
| 1 | [Production Readiness Report](PRODUCTION_READINESS_REPORT.md) | 25-section PRR covering capabilities, architecture, data inventory, test suite, security, and deployment checklists | Engineering, QA, Leadership |
| 2 | [Project Bible](PROJECT_BIBLE.md) | Single source of truth: architecture, directory structure, data models, dependencies, naming conventions, quality gates | Engineering |
| 3 | [Architecture Guide](ARCHITECTURE_GUIDE.md) | System architecture, clinical scale calculator design, stroke triage pipeline, ATN staging, multi-collection search, workflow engine | Engineering, Clinical Informatics |
| 4 | [White Paper](WHITE_PAPER.md) | Problem statement (3B affected, data fragmentation), RAG solution, clinical capabilities, results | External, Publications |
| 5 | [Deployment Guide](DEPLOYMENT_GUIDE.md) | Standalone Docker, integrated HCLS AI Factory, local development, configuration reference, troubleshooting | DevOps, Engineering |
| 6 | [Demo Guide](DEMO_GUIDE.md) | 5 clinical demo scenarios: acute stroke, memory clinic AD, drug-resistant epilepsy, brain mass, MS monitoring | Sales, Demo, GTC |
| 7 | [Learning Guide -- Foundations](LEARNING_GUIDE_FOUNDATIONS.md) | Neuroscience primer: brain anatomy, stroke types, dementia spectrum, seizure classification, MS pathology, clinical scales | New Team Members, Non-Clinical |
| 8 | [Learning Guide -- Advanced](LEARNING_GUIDE_ADVANCED.md) | ATN framework, DAWN/DEFUSE-3 criteria, NEDA-3, EMG/NCS patterns, DBS candidacy, epilepsy surgery, molecular neuro-oncology | Clinical Engineers, Domain Experts |
| 9 | [Research Paper](NEUROLOGY_INTELLIGENCE_AGENT_RESEARCH_PAPER.md) | Technical research paper on agent design and implementation | Academic, Publications |

---

## Quick Links

- **API Health:** `GET http://localhost:8528/health`
- **API Docs (Swagger):** `http://localhost:8528/docs`
- **Streamlit UI:** `http://localhost:8529`
- **GitHub Repository:** `ai_agent_adds/neurology_intelligence_agent/`

---

## Key Statistics

| Metric | Value |
|---|---|
| Milvus collections | 14 |
| Clinical workflows | 9 (8 domain-specific + general) |
| Clinical scale calculators | 10 |
| Disease domains | 10 |
| Drugs cataloged | 43 |
| Genes tracked | 38 |
| Entity aliases | 251+ |
| Synonym maps | 16 |
| Automated tests | 209 across 12 modules |
| Estimated vector records | 855,000 |

---

*HCLS AI Factory / GTC Europe 2026*
