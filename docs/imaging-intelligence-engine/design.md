# Clinical Imaging Engine -- Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 2.1.0
**License:** Apache 2.0

---

## 1. Purpose

This document describes the high-level design of the Clinical Imaging Engine (Engine 4), a multi-modal medical imaging analysis system that integrates 9 NVIDIA NIM clients, agentic reasoning (AIQ Plan/Execute/Reflect/Refine), cross-modal reasoning, NeMo Guardrails, and RAG-powered clinical interpretation. The engine delivers 9 clinical workflows, 13 Milvus collections (38,028 vectors including 1,938 real PubMed papers), 1,324 tests, and deploys on a single NVIDIA DGX Spark ($4,699) using 20 NVIDIA technologies (Community Edition, all free).

## 2. Design Goals

1. **Multi-modality support** -- CT, MRI, X-ray, ultrasound, mammography, and pathology imaging
2. **NIM integration** -- 9 NVIDIA NIM clients (VISTA-3D, MAISI, VILA-M3, LLM, NV-Segment-CT, Nemotron Nano, NV-Generate-CT, NV-Generate-MR, NV-Reason-CXR stub) with 3-tier fallback (cloud, local, mock)
3. **Agentic reasoning** -- AIQ Plan/Execute/Reflect/Refine with 6 tools
4. **Cross-modal analysis** -- Correlate imaging findings with genomic, clinical, and literature data via 8 cross-modal triggers
5. **Guardrails** -- NeMo Guardrails for PII protection, evidence grounding, and disclaimer enforcement
6. **Radiomics** -- ~1,500 radiomics features via PyRadiomics-CUDA
7. **Report NLP** -- Full radiology report parsing pipeline
8. **Streaming** -- Holoscan real-time ultrasound/endoscopy pipeline
9. **FLARE integration** -- Federated learning for medical image segmentation with MONAI Label interactive annotation
10. **Platform integration** -- Operates within the HCLS AI Factory ecosystem as Engine 4

## 3. Architecture Overview

- **API Layer** (FastAPI) -- Imaging analysis endpoints, cross-modal queries, report generation
- **Intelligence Layer** -- 9 NIM clients, AIQ agentic reasoning with 6 tools, RAG retrieval, cross-modal reasoning, NeMo Guardrails
- **Data Layer** (Milvus) -- 13 vector collections (38,028 vectors) for radiology literature, imaging protocols, guidelines, radiomics, and reports
- **Presentation Layer** (React Portal + Streamlit) -- Full React portal with 10 pages, interactive Streamlit workbench with NVIDIA dark theme, Three.js rotating point cloud visualization
- **Streaming Layer** (Holoscan) -- Real-time ultrasound/endoscopy pipeline
- **Deploy Layer** (MONAI Deploy) -- 9 MAPs packaged for clinical deployment

For detailed technical architecture, see [ARCHITECTURE_GUIDE.md](ARCHITECTURE_GUIDE.md).

## 4. Key Design Decisions

| Decision | Rationale |
|---|---|
| 3-tier NIM fallback | Cloud -> local -> mock ensures availability across deployment targets |
| 9 NIM clients | Comprehensive coverage: segmentation, generation, reasoning, language |
| AIQ agentic reasoning | Plan/Execute/Reflect/Refine with 6 tools for multi-step clinical analysis |
| NeMo Guardrails | PII protection, evidence grounding, disclaimer enforcement |
| Cross-modal correlation | 8 triggers bridging imaging findings with genomic and clinical context |
| FLARE federated learning | Privacy-preserving model training across institutions with MONAI Label |
| 3-tier deployment | Community/Enterprise/Research deployment models |
| Mock mode for development | Full API compatibility without GPU requirements |

## 5. Disclaimer

This system is a research and decision-support tool. It is not FDA-cleared or CE-marked and is not intended for independent clinical decision-making. All outputs should be reviewed by qualified clinical professionals.

---

*Clinical Imaging Engine -- Design Document v2.1.0*
*HCLS AI Factory -- Apache 2.0 | Author: Adam Jones | March 2026*
