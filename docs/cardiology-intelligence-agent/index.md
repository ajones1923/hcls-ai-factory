# Cardiology Intelligence Agent — Documentation Index

> **Part of the [Precision Intelligence Network](../engines/precision-intelligence.md)** — one of 11 specialized agents sharing a common molecular foundation within the HCLS AI Factory.

RAG-powered cardiovascular clinical decision support system built on Milvus, Claude, and BGE-small-en-v1.5.
Part of the [HCLS AI Factory](https://hcls-ai-factory.org) precision medicine platform.

**Version:** 2.0.0 | **Author:** Adam Jones | **Date:** March 2026

---

## Documentation

| Document | Description |
|----------|-------------|
| [Project Bible](PROJECT_BIBLE.md) | Complete reference — collections, workflows, API, configuration |
| [Architecture Guide](ARCHITECTURE_GUIDE.md) | System design, data flow, component interactions |
| [White Paper](WHITE_PAPER.md) | Clinical rationale, validation, DGX Spark performance |
| [Deployment Guide](DEPLOYMENT_GUIDE.md) | Docker Compose, manual setup, Milvus tuning, security |
| [Demo Guide](DEMO_GUIDE.md) | 10-tab walkthrough, sample queries, 5 demo scenarios |
| [Learning Guide — Foundations](LEARNING_GUIDE_FOUNDATIONS.md) | Cardiology primer for engineers and data scientists |
| [Learning Guide — Advanced](LEARNING_GUIDE_ADVANCED.md) | Cardiac MRI, hemodynamics, channelopathy genetics |
| [Production Readiness Report](PRODUCTION_READINESS_REPORT.md) | Capability matrix, test suite, deployment checklist |
| [Research Paper](CARDIOLOGY_INTELLIGENCE_AGENT_RESEARCH_PAPER.md) | Academic-style paper with product requirements |

## Quick Links

- **Streamlit UI:** Port 8536
- **FastAPI API:** Port 8126
- **11 Clinical Workflows** | **6 Risk Calculators** | **13 Vector Collections**
- **1,966 Tests** | 100% Pass Rate | <1.2s Runtime

## Platform

Built for NVIDIA DGX Spark as part of the HCLS AI Factory:
Patient DNA → Drug Candidates in <5 hours.

---

*Apache 2.0 License*

---

!!! warning "Clinical Decision Support Disclaimer"
    This agent is a clinical decision support research tool. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
