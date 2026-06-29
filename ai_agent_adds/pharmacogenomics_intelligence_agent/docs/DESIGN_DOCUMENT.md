# Pharmacogenomics Intelligence Agent -- Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.3.0
**License:** Apache 2.0

---

## 1. Purpose

This document describes the high-level design of the Pharmacogenomics Intelligence Agent, a RAG-powered clinical decision support system for pharmacogenomic variant interpretation, drug-gene interaction analysis, and genotype-guided prescribing.

## 2. Design Goals

1. **Genotype-guided prescribing** -- CPIC guideline-driven dosing recommendations
2. **Multi-source evidence** -- CPIC, PharmVar, PharmGKB, FDA labeling, PubMed integration
3. **VCF-to-recommendation pipeline** -- Direct variant interpretation with star allele mapping
4. **Evidence-grounded responses** -- All recommendations backed by CPIC level of evidence
5. **Platform integration** -- Operates within the HCLS AI Factory ecosystem

## 3. Architecture Overview

- **API Layer** (FastAPI, port 8107) -- Clinical endpoints, genotype interpretation, drug interaction queries
- **Intelligence Layer** -- Multi-collection RAG retrieval with pharmacogenomic-specific filtering
- **Data Layer** (Milvus) -- Vector collections for PGx literature, CPIC guidelines, drug labels, clinical trials
- **Presentation Layer** (Streamlit, port 8507) -- Interactive pharmacogenomics dashboard

For detailed technical architecture, see [ARCHITECTURE_GUIDE.md](ARCHITECTURE_GUIDE.md).

## 4. Key Design Decisions

| Decision | Rationale |
|---|---|
| CPIC guideline integration | Gold standard for pharmacogenomic clinical recommendations |
| Star allele resolution | PharmVar-based allele definitions for accurate genotype calls |
| cyvcf2 VCF parsing | High-performance C-backed VCF parsing for real-time analysis |
| Multi-source knowledge base | Combined CPIC, PharmGKB, FDA, and PubMed for comprehensive coverage |

## 5. Disclaimer

This system is a research and decision-support tool. It is not FDA-cleared or CE-marked and is not intended for independent clinical decision-making. All outputs should be reviewed by qualified clinical professionals.

---

*Pharmacogenomics Intelligence Agent -- Design Document v1.3.0*
*HCLS AI Factory -- Apache 2.0 | Author: Adam Jones | March 2026*
