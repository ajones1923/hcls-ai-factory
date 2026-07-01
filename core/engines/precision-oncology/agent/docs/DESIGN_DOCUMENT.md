# Precision Oncology Intelligence Agent -- Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.3.0
**License:** Apache 2.0

---

## 1. Purpose

This document describes the high-level design of the Precision Oncology Intelligence Agent, a closed-loop clinical decision support system that transforms variant call format (VCF) files into molecular tumor board (MTB) packets using RAG-powered evidence synthesis.

## 2. Design Goals

1. **End-to-end VCF-to-MTB workflow** -- Automated variant interpretation, therapy matching, and report generation
2. **Evidence-grounded responses** -- Every recommendation backed by ClinVar, OncoKB, NCCN, and PubMed citations
3. **Modular architecture** -- Pluggable data sources, LLM providers, and output formats
4. **Clinical safety** -- Confidence scoring, disclaimer generation, and audit trails
5. **Platform integration** -- Seamless operation within the HCLS AI Factory three-stage pipeline

## 3. Architecture Overview

The agent follows a four-layer architecture:

- **API Layer** (FastAPI) -- RESTful endpoints for queries, reports, and clinical workflows
- **Intelligence Layer** -- Multi-collection RAG retrieval, comparative analysis, variant interpretation
- **Data Layer** (Milvus) -- Vector collections for oncology literature, clinical trials, drug data, guidelines
- **Presentation Layer** (Streamlit) -- Interactive clinical dashboard with NVIDIA dark theme

For detailed technical architecture, see [ARCHITECTURE_GUIDE.md](ARCHITECTURE_GUIDE.md).

## 4. Key Design Decisions

| Decision | Rationale |
|---|---|
| Multi-collection Milvus | Separate collections per knowledge domain for filtered, domain-specific retrieval |
| Claude as LLM | Strong clinical reasoning, structured output, citation adherence |
| BGE-small-en-v1.5 embeddings | Balanced performance/size for biomedical text (384-dim) |
| FastAPI + Streamlit separation | Independent scaling of API and UI components |
| Pydantic v2 schemas | Strict validation of clinical data models |

## 5. Data Flow

```
VCF Input --> Variant Extraction --> Multi-Collection RAG Retrieval
    --> Evidence Synthesis (Claude) --> Confidence Scoring
    --> MTB Report Generation --> Output (Markdown/JSON/PDF/FHIR)
```

## 6. Integration Points

- **HCLS AI Factory Stage 1** -- Receives VCF from genomics pipeline
- **HCLS AI Factory Stage 3** -- Feeds druggable targets to drug discovery pipeline
- **External APIs** -- ClinicalTrials.gov, PubMed, OncoKB

## 7. Disclaimer

This system is a research and decision-support tool. It is not FDA-cleared or CE-marked and is not intended for independent clinical decision-making. All outputs should be reviewed by qualified clinical professionals.

---

*Precision Oncology Intelligence Agent -- Design Document v1.3.0*
*HCLS AI Factory -- Apache 2.0 | Author: Adam Jones | March 2026*
