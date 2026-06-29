# Neurology Intelligence Agent -- Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.3.0
**License:** Apache 2.0

---

## 1. Purpose

This document describes the high-level design of the Neurology Intelligence Agent, a RAG-powered neurology clinical decision support system providing 8 clinical workflows and 10 validated clinical scale calculators.

## 2. Design Goals

1. **Comprehensive neurology workflows** -- Acute stroke, dementia, epilepsy, brain tumors, MS, Parkinson's, headache, neuromuscular evaluation
2. **Validated scale calculators** -- NIHSS, GCS, MoCA, MDS-UPDRS Part III, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn-Yahr
3. **Guideline integration** -- AAN, AHA/ASA, ILAE, ICHD-3, WHO CNS 2021, NCCN, McDonald 2017, MDS criteria
4. **Evidence-grounded responses** -- All recommendations backed by citations
5. **Platform integration** -- Operates within the HCLS AI Factory ecosystem

## 3. Architecture Overview

- **API Layer** (FastAPI, port 8528) -- Clinical endpoints, scale calculators, report generation
- **Intelligence Layer** -- Multi-collection RAG retrieval with neurology-specific filtering
- **Data Layer** (Milvus) -- Vector collections for neurology literature, guidelines, trials
- **Presentation Layer** (Streamlit, port 8529) -- Interactive neurology dashboard

For detailed technical architecture, see [ARCHITECTURE_GUIDE.md](ARCHITECTURE_GUIDE.md).

## 4. Key Design Decisions

| Decision | Rationale |
|---|---|
| 8 specialized workflows | Domain-specific retrieval strategies per neurological condition |
| 10 validated scale calculators | Peer-reviewed scoring algorithms with clinical validation |
| Pydantic settings with NEURO_ prefix | Clean environment variable namespacing |
| SSE event streaming | Real-time progressive output for clinical workflows |

## 5. Disclaimer

This system is a research and decision-support tool. It is not FDA-cleared or CE-marked and is not intended for independent clinical decision-making. All outputs should be reviewed by qualified clinical professionals.

---

*Neurology Intelligence Agent -- Design Document v1.3.0*
*HCLS AI Factory -- Apache 2.0 | Author: Adam Jones | March 2026*
