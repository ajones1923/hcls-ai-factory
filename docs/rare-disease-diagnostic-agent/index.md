# Rare Disease Diagnostic Agent -- Documentation Index


![Rare Disease Diagnostic Agent Architecture Infographic](infographic.jpg)

*Source: [github.com/ajones1923/rare-disease-diagnostic-agent](https://github.com/ajones1923/rare-disease-diagnostic-agent)*

> **Part of the [Precision Intelligence Network](../engines/precision-intelligence.md)** — one of 11 specialized agents sharing a common molecular foundation within the HCLS AI Factory.

RAG-powered rare disease diagnostic support system built on Milvus, Claude, and BGE-small-en-v1.5.
Part of the [HCLS AI Factory](https://hcls-ai-factory.org) precision medicine platform.

**Author:** Adam Jones | **Date:** March 2026

---

## Documentation

| Document | Description |
|----------|-------------|
| [Project Bible](PROJECT_BIBLE.md) | Complete reference -- collections, workflows, API, configuration |
| [Architecture Guide](ARCHITECTURE_GUIDE.md) | System design, data flow, HPO matching pipeline, ACMG pipeline |
| [White Paper](WHITE_PAPER.md) | Diagnostic odyssey problem, RAG-based solution, validation |
| [Deployment Guide](DEPLOYMENT_GUIDE.md) | Docker Compose, manual setup, Milvus tuning, security |
| [Demo Guide](DEMO_GUIDE.md) | 5-tab walkthrough, sample HPO queries, ACMG demo, gene therapy demo |
| [Learning Guide -- Foundations](LEARNING_GUIDE_FOUNDATIONS.md) | Rare disease primer, HPO ontology, ACMG criteria, inheritance, gene therapy |
| [Learning Guide -- Advanced](LEARNING_GUIDE_ADVANCED.md) | ACMG 28-criteria deep dive, phenotype matching, family segregation, natural history |
| [Production Readiness Report](PRODUCTION_READINESS_REPORT.md) | Capability matrix, test suite, deployment checklist |
| [Research Paper](RARE_DISEASE_DIAGNOSTIC_AGENT_RESEARCH_PAPER.md) | Full technical paper with PRD elements |

## Quick Links

- **Streamlit UI:** Port 8544
- **FastAPI API:** Port 8134
- **10 Diagnostic Workflows** | **14 Vector Collections** | **13 Disease Categories**
- **193 Tests** | 100% Pass Rate | <0.2s Runtime

---

*Apache 2.0 License*

---

!!! warning "Clinical Decision Support Disclaimer"
    This agent is a clinical decision support research tool. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
