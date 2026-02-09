---
title: "Stage 2: Evidence RAG"
search:
  boost: 2
tags:
  - RAG
  - Milvus
  - ClinVar
  - AlphaMissense
  - Claude
  - Vector Database
---

# Stage 2: Evidence RAG

<div class="stage-hero stage-hero-teal">
  <div class="stage-number">02</div>
  <div class="stage-intro">
    <p class="stage-tagline">Finding the needle in 11.7 million variants</p>
    <p class="stage-time">Interactive queries in seconds</p>
  </div>
</div>

---

## What This Stage Does

Stage 1 produced **11.7 million variants**. Most are harmless — they're just what makes you, you. The challenge is finding the ones that actually **cause disease** and can be **targeted with drugs**.

Stage 2 uses **Retrieval-Augmented Generation (RAG)** to:

1. **Annotate** — Cross-reference every variant against clinical databases
2. **Embed** — Convert variants into searchable vector representations
3. **Reason** — Use AI to interpret evidence and identify targets

---

## The Funnel

This is the heart of what the platform does — narrowing millions of variants to actionable targets:

```mermaid
flowchart TD
    A["11,700,000 variants discovered"]
    B["3,500,000 pass quality filters"]
    C["35,616 match ClinVar records"]
    D["6,831 have pathogenicity predictions"]
    E["2,400 high-impact, disease-causing"]
    F["847 in druggable genes"]

    A -->|"Quality filter (QUAL>30)"| B
    B -->|"ClinVar annotation"| C
    C -->|"AlphaMissense scoring"| D
    D -->|"Impact classification"| E
    E -->|"Druggability filter"| F

    style A fill:#B3E5FC,stroke:#0288D1
    style B fill:#81D4FA,stroke:#0288D1
    style C fill:#4FC3F7,stroke:#0277BD
    style D fill:#29B6F6,stroke:#0277BD,color:#fff
    style E fill:#039BE5,stroke:#01579B,color:#fff
    style F fill:#00B4D8,stroke:#0077B6,color:#fff
```

From 11.7 million down to 847. That's the power of AI-driven annotation.

---

## Knowledge Base

The RAG system is grounded in curated biomedical knowledge:

| Source | Coverage |
|--------|----------|
| **ClinVar** | 4.1 million clinically studied variants |
| **AlphaMissense** | DeepMind's pathogenicity predictions |
| **Gene Panel** | 201 genes across 13 therapeutic areas |
| **Druggability** | 171 genes (85%) with known drug targets |

### Therapeutic Areas Covered

Neurology · Oncology · Cardiovascular · Metabolic · Immunology · Rare Disease · Ophthalmology · Dermatology · Respiratory · Hematology · Musculoskeletal · Endocrine · Infectious Disease

---

## How RAG Works

**Retrieval-Augmented Generation** combines search with AI reasoning:

```mermaid
flowchart LR
    Q["Natural Language\nQuery"]
    E["BGE Embedding\n384 dimensions"]
    M["Milvus Vector Search\n3.5M variants"]
    R["Retrieved Evidence\nClinVar + AlphaMissense"]
    C["Claude AI\nGrounded Reasoning"]
    A["Grounded Answer\nwith citations"]

    Q --> E --> M --> R --> C --> A

    style Q fill:#E1BEE7,stroke:#7B1FA2
    style E fill:#B3E5FC,stroke:#0288D1
    style M fill:#00B4D8,stroke:#0077B6,color:#fff
    style R fill:#B2DFDB,stroke:#00796B
    style C fill:#FFE082,stroke:#FFA000
    style A fill:#C8E6C9,stroke:#388E3C
```

1. **You ask a question** in natural language
2. **Vector search** finds the most relevant variants and annotations
3. **Evidence is retrieved** from the knowledge base
4. **Claude reasons** over the evidence to provide grounded answers

The AI can only cite variants and evidence that **actually exist in this patient's data**. It's not hallucinating — it's reasoning over real genomic evidence.

---

## Example Query

**Question:** "What pathogenic variants does this patient have in genes associated with neurodegeneration?"

**Response:** The system identifies a VCP variant (chr9:35,065,263 G>A) with:

- ClinVar classification: **Pathogenic**
- AlphaMissense score: **0.87** (threshold: 0.564)
- Associated disease: **Frontotemporal Dementia**
- Druggability: **Yes** — VCP is a known therapeutic target

This variant becomes the input for [Stage 3: Drug Discovery](stage-3-drug-discovery.md).

---

## Technology Stack

<div class="tech-pills">

- **Milvus 2.4** — Vector database for 3.5M variant embeddings
- **Anthropic Claude** — AI reasoning over retrieved evidence
- **ClinVar** — NIH clinical variant database
- **AlphaMissense** — DeepMind pathogenicity predictions
- **Streamlit** — Interactive chat interface
- **Flask** — REST API for programmatic access

</div>

---

## By the Numbers

| Metric | Value |
|--------|-------|
| Variants indexed | 3.5 million |
| Embedding dimensions | 384 |
| Query latency | < 2 seconds |
| Knowledge sources | 4 (ClinVar, AlphaMissense, Gene Panel, PDB) |
| Genes covered | 201 |
| Druggable targets | 171 (85%) |

---

## Why This Matters

Traditional variant interpretation requires:

- **Manual literature review** — Hours per variant
- **Expert curation** — Scarce clinical geneticists
- **Fragmented tools** — Separate databases, no unified view

The RAG pipeline delivers:

- **Instant answers** — Seconds, not hours
- **Grounded reasoning** — Every claim traceable to evidence
- **Unified view** — All annotations in one place
- **Natural language** — No bioinformatics expertise required

