---
title: Precision Intelligence Network
description: 11 domain-specialized AI agents sharing a common molecular foundation
---

# Precision Intelligence Network

**Stage 2 of the HCLS AI Factory Pipeline**

The Precision Intelligence Network is a coordinated system of 11 domain-specialized intelligence agents, each built on RAG (Retrieval-Augmented Generation) architecture with Milvus vector search and Claude AI. All 11 agents share read-only access to the `genomic_evidence` collection (3.56M vectors) produced by the [Genomic Foundation Engine](genomic-foundation.md), while maintaining their own domain-specific vector collections.

## What It Does

Raw genomic variants from Stage 1 become clinically actionable intelligence through domain-specific interpretation, cross-agent coordination, and evidence-grounded reasoning. Each agent transforms the shared molecular foundation into specialized insights for its clinical domain.

## The 11 Intelligence Agents

| Agent | Domain | UI Port | API Port | Collections | Key Capabilities |
|-------|--------|---------|----------|-------------|------------------|
| [Precision Oncology](../precision-oncology-agent/index.md) | Cancer | 8526 | 8527 | 11 | MTB packets, therapy ranking, trial matching, FHIR R4 |
| [CAR-T Intelligence](../cart-intelligence-agent/index.md) | Cell Therapy | 8521 | -- | 11 | Cross-collection evidence, comparative analysis, deep research |
| [Imaging Intelligence](../imaging-intelligence-agent/index.md) | Medical Imaging | 8525 | 8524 | 10 | VISTA-3D, MAISI, VILA-M3 workflows, FHIR R4 export |
| [Precision Biomarker](../precision-biomarker-agent/index.md) | Biomarkers | 8528 | 8529 | 14 | PhenoAge/GrimAge, 9-domain risk, genotype-aware interpretation |
| [Pharmacogenomics](../pharmacogenomics-intelligence-agent/index.md) | Drug-Gene | 8507 | 8107 | 15 | Star allele calling, CPIC guidelines, 9 dosing algorithms |
| [Precision Autoimmune](../precision-autoimmune-agent/index.md) | Autoimmune | 8531 | 8532 | 14 | Autoantibody interpretation, HLA analysis, flare prediction |
| [Neurology Intelligence](../neurology-intelligence-agent/index.md) | Neurology | 8529 | 8528 | 14 | Stroke triage, dementia evaluation, EDSS scoring |
| [Cardiology Intelligence](../cardiology-intelligence-agent/index.md) | Cardiovascular | 8536 | 8126 | 13 | 11 clinical workflows, 6 risk calculators |
| [Clinical Trial Intelligence](../clinical-trial-intelligence-agent/index.md) | Clinical Trials | 8128 | 8538 | 14 | Trial optimization, adaptive design, biomarker strategy |
| [Rare Disease Diagnostic](../rare-disease-diagnostic-agent/index.md) | Rare Disease | 8544 | 8134 | 14 | HPO matching, ACMG classification, gene therapy tracking |
| [Single-Cell Intelligence](../single-cell-intelligence-agent/index.md) | Single-Cell | 8130 | 8540 | 12 | Cell type annotation, TME profiling, drug response |

## Shared Molecular Foundation

All 11 agents connect to the `genomic_evidence` collection as a read-only data source. This shared foundation ensures:

- **Consistency** -- Every agent reasons from the same molecular data
- **Efficiency** -- 3.56M vectors are indexed once, queried by all
- **Traceability** -- Variant-level lineage from FASTQ through every agent's output

## Cross-Agent Coordination

Agents communicate through Server-Sent Events (SSE) and cross-modal triggers:

- **Imaging → Genomics**: Lung-RADS 4A+ triggers EGFR/ALK/ROS1/KRAS genomic queries
- **Oncology → Drug Discovery**: Actionable variants trigger molecule generation in Stage 3
- **Biomarker → Pharmacogenomics**: Genotype findings trigger drug interaction checks
- **Autoimmune → Biomarker**: HLA associations trigger biomarker monitoring recommendations

## Shared Technology Stack

| Component | Technology |
|-----------|-----------|
| Vector Database | Milvus 2.4 with IVF_FLAT/COSINE indexes |
| Embeddings | BGE-small-en-v1.5 (384-dim) |
| LLM | Claude (Anthropic API) |
| UI Framework | Streamlit with NVIDIA dark theme |
| API Framework | FastAPI |
| Export Formats | Markdown, JSON, PDF, FHIR R4 |
| Hardware Target | NVIDIA DGX Spark |

## Agent Architecture Pattern

Every agent follows a consistent architecture:

```
User Query / Patient Data
    |
    v
[Domain-Specific Clinical Engines]
(deterministic scoring, classification, risk stratification)
    |
    v
[Multi-Collection RAG Search]
(BGE-small embedding → parallel Milvus search across N collections)
    |
    v
[Knowledge Graph Augmentation]
(domain entities, aliases, expansion maps)
    |
    v
[Claude LLM Synthesis]
(grounded response with citations)
    |
    v
[Clinical Output]
(Markdown | JSON | PDF | FHIR R4)
```

## Key Numbers

| Metric | Value |
|--------|-------|
| Total agents | 11 |
| Total domain-specific collections | 142 |
| Shared genomic vectors | 3.56M |
| Total test cases | 6,000+ |
| Therapeutic areas covered | 13 |
| Genes covered | 201 |
| Druggable targets | 171 |
| Query response time | <5 seconds (vector search) |

## Getting Started

- [Deployment Guide](../HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md) -- Full installation and configuration
- [Architecture](../architecture.md) -- System architecture overview
- [RAG/Chat Pipeline](../rag-chat-pipeline/README.md) -- Stage 2 pipeline details

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Intelligence Network agents are clinical decision support research tools. They are not FDA-cleared and are not intended as standalone diagnostic devices. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
