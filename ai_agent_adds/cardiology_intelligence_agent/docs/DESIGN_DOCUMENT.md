# Cardiology Intelligence Agent -- Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.3.0
**License:** Apache 2.0

---

## 1. Purpose

This document describes the high-level design of the Cardiology Intelligence Agent, a RAG-powered cardiovascular clinical decision support system that synthesizes cardiac imaging, electrophysiology, hemodynamics, and guideline evidence into actionable clinical recommendations.

## 2. Design Goals

1. **Guideline-aligned recommendations** -- ACC/AHA/ESC evidence-based clinical decision support
2. **Multi-domain cardiovascular intelligence** -- Heart failure, valvular disease, electrophysiology, interventional, preventive cardiology, cardio-oncology
3. **Validated clinical calculators** -- HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE II
4. **Evidence-grounded responses** -- Citations from ACC/AHA, ESC, PubMed, and ClinicalTrials.gov
5. **Platform integration** -- Operates within the HCLS AI Factory ecosystem

## 3. Architecture Overview

- **API Layer** (FastAPI) -- Clinical endpoints, scale calculators, report generation
- **Intelligence Layer** -- Multi-collection RAG retrieval with guideline-specific filtering
- **Data Layer** (Milvus) -- Vector collections for cardiovascular literature, guidelines, trials, drug data
- **Presentation Layer** (Streamlit) -- Interactive cardiovascular dashboard

For detailed technical architecture, see [ARCHITECTURE_GUIDE.md](ARCHITECTURE_GUIDE.md).

## 4. Key Design Decisions

| Decision | Rationale |
|---|---|
| Dedicated clinical scale endpoints | Validated scoring algorithms separate from LLM inference |
| Multi-guideline retrieval | Simultaneous ACC/AHA and ESC guideline search for comprehensive coverage |
| BGE-small-en-v1.5 embeddings | Optimized for biomedical text at 384 dimensions |
| Streaming SSE responses | Real-time progressive output for clinical workflows |

## 5. Disclaimer

This system is a research and decision-support tool. It is not FDA-cleared or CE-marked and is not intended for independent clinical decision-making. All outputs should be reviewed by qualified clinical professionals.

---

*Cardiology Intelligence Agent -- Design Document v1.3.0*
*HCLS AI Factory -- Apache 2.0 | Author: Adam Jones | March 2026*
