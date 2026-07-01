# Single-Cell Intelligence Agent -- Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.3.0
**License:** Apache 2.0

---

## 1. Purpose

This document describes the high-level design of the Single-Cell Intelligence Agent, a RAG-powered system for single-cell RNA sequencing analysis, cell type annotation, spatial transcriptomics, tumor microenvironment characterization, and CAR-T target validation.

## 2. Design Goals

1. **Cell type annotation** -- Automated classification across 32 cell types with marker gene validation
2. **Spatial transcriptomics** -- Support for 4 spatial platforms (10x Visium, MERFISH, Slide-seq, Xenium)
3. **TME characterization** -- Tumor microenvironment classification and immune cell profiling
4. **CAR-T target validation** -- Target expression analysis for cell therapy development
5. **Platform integration** -- Operates within the HCLS AI Factory ecosystem

## 3. Architecture Overview

- **API Layer** (FastAPI, port 8130) -- Clinical endpoints, workflow dispatch, knowledge queries
- **Intelligence Layer** -- Multi-collection RAG retrieval, cell type classification, spatial analysis
- **Data Layer** (Milvus) -- Vector collections for single-cell literature, marker genes, drug data
- **Presentation Layer** (Streamlit, port 8540) -- Interactive single-cell analysis dashboard

For detailed technical architecture, see [ARCHITECTURE_GUIDE.md](ARCHITECTURE_GUIDE.md).

## 4. Key Design Decisions

| Decision | Rationale |
|---|---|
| 32 cell types, 55 markers | Curated reference set covering major tissue and tumor types |
| Multi-workflow dispatch | Generic `/workflow/{type}` endpoint for extensible analysis pipelines |
| 4 spatial platforms | Coverage of major spatial transcriptomics technologies |
| Knowledge versioning | Explicit version tracking for reproducible analyses |

## 5. Disclaimer

This system is a research and decision-support tool. It is not FDA-cleared or CE-marked and is not intended for independent clinical decision-making. All outputs should be reviewed by qualified clinical professionals.

---

*Single-Cell Intelligence Agent -- Design Document v1.3.0*
*HCLS AI Factory -- Apache 2.0 | Author: Adam Jones | March 2026*
