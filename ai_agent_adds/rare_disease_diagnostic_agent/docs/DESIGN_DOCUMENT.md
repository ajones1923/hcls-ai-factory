# Rare Disease Diagnostic Agent -- Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.3.0
**License:** Apache 2.0

---

## 1. Purpose

This document describes the high-level design of the Rare Disease Diagnostic Agent, a RAG-powered diagnostic support system for rare disease identification, phenotype-to-genotype mapping, and diagnostic odyssey reduction.

## 2. Design Goals

1. **Diagnostic support** -- Phenotype-driven differential diagnosis for rare diseases
2. **Multi-format reporting** -- Markdown, JSON, PDF, FHIR R4, and GA4GH Phenopacket v2 output
3. **Knowledge integration** -- Orphanet, OMIM, HPO, ClinVar, and literature-backed evidence
4. **Phenopacket interoperability** -- GA4GH Phenopacket v2 for rare disease data exchange
5. **Platform integration** -- Operates within the HCLS AI Factory ecosystem

## 3. Architecture Overview

- **API Layer** (FastAPI, port 8134) -- Diagnostic endpoints, report generation, phenotype queries
- **Intelligence Layer** -- Multi-collection RAG retrieval, phenotype-genotype correlation
- **Data Layer** (Milvus) -- Vector collections for rare disease literature, phenotype ontologies, genetic data
- **Presentation Layer** (Streamlit, port 8544) -- Interactive diagnostic dashboard

For detailed technical architecture, see [ARCHITECTURE_GUIDE.md](ARCHITECTURE_GUIDE.md).

## 4. Key Design Decisions

| Decision | Rationale |
|---|---|
| GA4GH Phenopacket v2 | Standard format for rare disease phenotypic data exchange |
| FHIR R4 DiagnosticReport | HL7 interoperability for clinical system integration |
| HPO-based phenotype matching | Human Phenotype Ontology as standard phenotype vocabulary |
| Multi-format output | Clinical, research, and interoperability use cases |

## 5. Disclaimer

This system is a research and decision-support tool. It is not FDA-cleared or CE-marked and is not intended for independent clinical decision-making. All outputs should be reviewed by qualified clinical professionals.

---

*Rare Disease Diagnostic Agent -- Design Document v1.3.0*
*HCLS AI Factory -- Apache 2.0 | Author: Adam Jones | March 2026*
