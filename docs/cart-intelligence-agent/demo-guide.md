---
search:
  boost: 2
tags:
  - Demo
  - Walkthrough
  - CAR-T
  - Cell Therapy
  - Immunotherapy
---

# CAR-T Intelligence Agent — Demo Guide

> **Step-by-step walkthrough for demonstrating the CAR-T Intelligence Agent on DGX Spark.**
>
> License: Apache 2.0 | Date: February 2026

---

## Demo Overview

| Parameter | Value |
|---|---|
| **Route A Duration** | 15 minutes (standalone agent) |
| **Route B Duration** | 25-30 minutes (cross-platform integration) |
| **Hardware** | NVIDIA DGX Spark (GB10, 128 GB unified) |
| **Knowledge Base** | 6,266+ vectors across 10 owned collections + 3.56M genomic vectors |
| **Knowledge Graph** | 25 target antigens, 6 FDA products, 8 toxicity profiles |
| **LLM** | Claude Sonnet 4.6 (Anthropic) |
| **Export Formats** | Markdown, JSON, PDF |

### What the Audience Will See

1. A health dashboard showing 11 vector collections spanning the entire CAR-T development lifecycle
2. Cross-functional RAG queries pulling evidence from literature, trials, constructs, assays, and manufacturing — simultaneously
3. Automatic comparative analysis triggered by natural language — "Compare X vs Y" produces structured tables
4. 12-16 millisecond retrieval across 11 collections with 1,496-term query expansion
5. Knowledge graph augmentation with 25 target antigens, 6 FDA-approved products, and 8 toxicity profiles
6. Deep research mode: autonomous agent decomposing complex questions into sub-queries
7. Clickable PubMed and ClinicalTrials.gov citations grounding every answer
8. Professional reports exported as Markdown, JSON, and NVIDIA-themed PDF
9. *(Route B)* Cross-pipeline data sharing — CAR-T intelligence layered on 3.56 million genomic variant vectors
10. *(Route B)* The complete loop: Patient DNA → Variant Analysis → CAR-T Intelligence → Drug Candidates

---

## Pre-Demo Setup

### Step 1: Verify Hardware

```bash
# Verify DGX Spark GPU
nvidia-smi
# Expected: GB10 GPU, 128 GB unified memory

# Verify ARM64 architecture
uname -m
# Expected: aarch64
```

### Step 2: Set Environment Variables

```bash
cp .env.example .env

# Required variables:
# ANTHROPIC_API_KEY=sk-ant-...     (for Claude RAG synthesis)
```

### Step 3: Start Services

```bash
cd ai_agent_adds/cart_intelligence_agent
docker compose up -d

# This starts:
# Milvus (etcd + MinIO + standalone)
# Streamlit UI (port 8521)
# FastAPI server (port 8522)
```

### Step 4: Verify All Services Healthy

```bash
curl -s http://localhost:8522/health | python3 -m json.tool
```

Expected response:

```json
{
  "status": "healthy",
  "collections": {
    "cart_literature": 5047,
    "cart_trials": 973,
    "cart_constructs": 6,
    "cart_assays": 45,
    "cart_manufacturing": 30,
    "cart_safety": 40,
    "cart_biomarkers": 43,
    "cart_regulatory": 25,
    "cart_sequences": 27,
    "cart_realworld": 30,
    "genomic_evidence": 3561170
  },
  "total_vectors": 3567436
}
```

All 11 collections should have non-zero counts.

### Step 5: Verify Knowledge Graph

```bash
curl -s http://localhost:8522/knowledge/stats | python3 -m json.tool
```

Expected: 25 target antigens, 6 FDA-approved products, 8 toxicity profiles, 10 manufacturing processes, 15+ biomarkers.

---

## Route A: Standalone Agent Demo (15 minutes)

### Opening (1 minute)

**Talking points:**

- "This is the CAR-T Intelligence Agent — it breaks down data silos across the entire CAR-T cell therapy development lifecycle."
- "11 Milvus collections covering literature, clinical trials, FDA-approved constructs, assays, manufacturing, safety, biomarkers, regulatory, sequences, real-world evidence, and genomic variants."
- "Every query searches all data sources simultaneously and Claude synthesizes cross-functional insights with citations."

**Show:** Health endpoint — highlight 11 collections, 6,266+ owned vectors, 3.56M genomic vectors.

---

### RAG Knowledge Query (3 minutes)

```bash
curl -s -X POST http://localhost:8522/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Why do CD19 CAR-T therapies fail in relapsed B-ALL?",
    "target_gene": "CD19"
  }' | python3 -m json.tool
```

Expected response highlights:

```json
{
  "answer": "CD19 CAR-T therapy failure in relapsed B-ALL occurs through several mechanisms...",
  "sources": [
    {"collection": "cart_literature", "score": 0.89, "text_snippet": "Antigen loss observed in 28% of relapses..."},
    {"collection": "cart_assays", "score": 0.85, "text_snippet": "CD19-negative relapses via trogocytosis..."},
    {"collection": "cart_trials", "score": 0.82, "text_snippet": "ELIANA trial: 12-month EFS 73%..."},
    {"collection": "cart_constructs", "score": 0.80, "text_snippet": "Kymriah tisagenlecleucel 4-1BB..."}
  ],
  "follow_up_questions": [
    "What strategies address CD19 antigen loss?",
    "How do dual-targeting constructs reduce relapse rates?",
    "What biomarkers predict CD19 CAR-T failure?"
  ],
  "confidence": 0.87,
  "processing_time_ms": 24100
}
```

**Talking points:**

- "One question, four collections hit simultaneously — literature, assays, trials, and constructs."
- "The answer identifies specific failure mechanisms: antigen loss in 28% of relapses, lineage switch in 10% with KMT2A rearrangements, trogocytosis, and T-cell exhaustion."
- "Every claim is backed by a citation with cosine similarity scores. The top hit scored 0.89 — very high relevance."
- "Claude suggests three follow-up questions, each targeting a different aspect of the problem."

---

### Comparative Analysis (4 minutes)

```bash
curl -s -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Compare 4-1BB vs CD28 costimulatory domains for DLBCL"
  }' | python3 -m json.tool
```

**Talking points:**

- "Watch what happens when I say 'compare'. The engine automatically detects this is a comparative query."
- "It parses two entities — 4-1BB and CD28 — and resolves each against the knowledge graph."
- "Dual retrieval runs: entity A and entity B are searched separately, then results are merged into a structured prompt."
- "Claude produces a comparison table with advantages, limitations, and clinical context."

**Expected output structure:**

```
## Comparison: 4-1BB vs CD28

| Dimension | 4-1BB (CD137) | CD28 |
|---|---|---|
| Peak expansion | Slower (4-8 weeks) | Faster (1-2 weeks) |
| T-cell persistence | Extended (months-years) | Shorter (weeks-months) |
| Exhaustion markers | Lower PD-1/LAG-3/TIM-3 | Higher exhaustion |
| CRS incidence | Lower grade 3-4 | Higher grade 3-4 |
| FDA products | Kymriah, Tecartus, Breyanzi, Abecma, Carvykti | Yescarta |
| Clinical data | ELIANA, TRANSCEND, KarMMa | ZUMA-1 |

### Advantages
**4-1BB:** Sustained persistence, lower exhaustion, broader product portfolio
**CD28:** Faster clinical response, stronger initial expansion
```

- "This comparison was auto-generated from dual retrieval across 11 collections. The engine found evidence for both domains across literature, trials, constructs, and assays."
- "Five of six FDA-approved CAR-T products use 4-1BB. Only Yescarta uses CD28. That's a clear trend."

**Comparison types you can demo:**

| Query | Entities Resolved |
|---|---|
| "Compare CD19 vs BCMA" | Target antigen vs target antigen |
| "Kymriah versus Carvykti" | Product → CD19 vs BCMA |
| "Compare CRS and ICANS" | Toxicity vs toxicity |
| "4-1BB vs CD28 for persistence" | Costimulatory domain comparison |
| "Lentiviral vs transposon vectors" | Manufacturing comparison |

---

### Evidence-Only Search (2 minutes)

```bash
curl -s -X POST http://localhost:8522/search \
  -H "Content-Type: application/json" \
  -d '{
    "question": "BCMA CAR-T resistance mechanisms in multiple myeloma",
    "target_antigen": "BCMA",
    "top_k": 5
  }' | python3 -m json.tool
```

**Talking points:**

- "This is raw retrieval — no LLM synthesis. Just the vector search results."
- "11 collections searched in parallel. Total retrieval time: 12-16 milliseconds."
- "Each hit shows collection, relevance score, and text snippet. This is what grounds the LLM."
- "Query expansion kicked in: 'BCMA' expanded to include 'B-cell maturation antigen', 'RRMM', 'plasma cell neoplasm', and related terms."
- "The top hit from cart_assays describes biallelic BCMA loss in 29% of relapses — that's sBCMA decoy signaling."

**Key performance metric:**

| Operation | Latency |
|---|---|
| Embedding | < 5 ms |
| 11-collection parallel search | 12-16 ms |
| Query expansion + re-search | 8-12 ms |
| Merge + deduplicate + rank | < 1 ms |
| **Total retrieval** | **~25 ms** |

---

### Knowledge Graph Stats (1 minute)

```bash
curl -s http://localhost:8522/knowledge/stats | python3 -m json.tool
```

**Talking points:**

- "The knowledge graph contains 25 CAR-T target antigens — from well-established CD19 and BCMA to emerging targets like GPRC5D and Claudin18.2."
- "All 6 FDA-approved products are modeled with their complete specifications: target antigen, scFv origin, costimulatory domain, vector type, and approval timeline."
- "8 toxicity profiles cover CRS, ICANS, B-cell aplasia, HLH/MAS, cytopenias, TLS, GvHD, and on-target/off-tumor toxicity."
- "This structured knowledge is automatically injected into every query prompt — Claude knows the entire landscape."

---

### Deep Research Mode (2 minutes)

**Show:** Streamlit UI at http://localhost:8521

**Toggle:** Enable "Deep Research Mode" in the sidebar.

**Type this query:**
> "What are the key challenges in developing solid tumor CAR-T therapies and what strategies are being explored?"

**Talking points:**

- "Deep Research Mode activates the full agent reasoning pipeline."
- "Step 1: Plan — the agent identifies this as a broad, multi-faceted question requiring decomposition."
- "Step 2: Decompose — it breaks the question into sub-queries: tumor microenvironment, antigen heterogeneity, trafficking, exhaustion, safety."
- "Step 3: Search — each sub-query gets its own retrieval pass across 11 collections."
- "Step 4: Evaluate — the agent checks evidence quality. If insufficient, it runs supplementary searches."
- "Step 5: Synthesize — Claude produces a comprehensive, multi-section answer with citations from every sub-query."
- "This is autonomous reasoning — the agent decides what to search for and how to structure the answer."

---

### Report Export (2 minutes)

```bash
# Markdown report
curl -s -X POST http://localhost:8522/api/report/markdown \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What manufacturing parameters predict CAR-T clinical response?",
    "target_gene": "CD19"
  }' | head -50

# JSON report
curl -s -X POST http://localhost:8522/api/report/json \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What manufacturing parameters predict CAR-T clinical response?"
  }' | python3 -m json.tool | head -30

# PDF report
curl -s -X POST http://localhost:8522/api/report/pdf \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What manufacturing parameters predict CAR-T clinical response?"
  }' --output cart_manufacturing_report.pdf
```

**Talking points:**

- "Three export formats: Markdown for sharing, JSON for programmatic consumption, PDF for clinical documentation."
- "The PDF is NVIDIA-themed with professional styling via ReportLab."
- "Every report includes: query metadata, RAG-synthesized analysis, evidence sources by collection, knowledge graph context, and citation links."
- "Download the PDF — open it and you'll see a publication-ready report with clickable PubMed links."

---

### Closing Route A (1 minute)

**Talking points:**

- "11 collections, 6,266 curated vectors, 25 target antigens, 12-16ms retrieval, automatic comparative analysis."
- "This is a complete CAR-T cell therapy intelligence platform — from target discovery through manufacturing to clinical safety."
- "Every answer is grounded in published evidence with clickable citations. No hallucination."

---

## Route B: Cross-Platform Integration Demo (25-30 minutes)

> **Prerequisite:** Complete Route A first, or start from a fresh session with all services running.
>
> This route demonstrates how the CAR-T Intelligence Agent connects to the full HCLS AI Factory — bridging cell therapy intelligence with genomic variant analysis and AI-driven drug discovery.

### Platform Overview (2 minutes)

**Show:** HCLS AI Factory landing page at http://localhost:8080

**Talking points:**

- "The HCLS AI Factory is a 3-stage precision medicine platform: GPU-accelerated genomics, RAG-grounded target identification, and AI-driven drug discovery."
- "The CAR-T Intelligence Agent extends this platform with specialized cell therapy knowledge."
- "All agents share the same Milvus vector database, BGE embeddings, and Claude LLM — creating a unified intelligence layer."

**Architecture overview:**

```
Stage 1: GPU Genomics (Parabricks)
    ↓ 11.7M variants → 3.56M quality-filtered
Stage 2: RAG Target Identification (Claude + Milvus)
    ↓ genomic_evidence: 3.56M vectors (shared)
CAR-T Intelligence Agent (11 collections)
    ↓ Target validation + construct intelligence
Stage 3: Drug Discovery (BioNeMo)
    ↓ 100 ranked drug candidates
```

---

### Shared Genomic Data Layer (3 minutes)

```bash
# Show the CAR-T agent's view of all collections
curl -s http://localhost:8522/collections | python3 -m json.tool
```

**Talking points:**

- "The CAR-T agent has 10 specialized collections — 6,266 curated vectors covering literature, trials, constructs, assays, manufacturing, safety, biomarkers, regulatory, sequences, and real-world evidence."
- "It also reads from `genomic_evidence` — 3,561,170 vectors from the genomics pipeline. These are the same variant annotations Stage 2 uses."
- "This is the key integration point. The same patient variants that Stage 2 analyzes for drug targets are available to inform CAR-T therapy decisions."

---

### Patient Variant Context (4 minutes)

```bash
curl -s -X POST http://localhost:8522/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What genomic variants in this patient affect CAR-T therapy response and toxicity prediction?"
  }' | python3 -m json.tool
```

**Talking points:**

- "This query bridges genomic variants and CAR-T intelligence. Claude searches both the therapy collections and the genomic evidence."
- "It identifies variants relevant to CRS prediction — ferritin, IL-6 pathway polymorphisms."
- "It flags variants that could affect T-cell expansion and persistence — immune checkpoint gene variants."
- "This is precision cell therapy — matching the patient's genome to the optimal CAR-T strategy."

---

### CAR-T + Genomics Bridge (4 minutes)

**Show:** Streamlit chat at http://localhost:8521

**Type this query:**
> "If the genomics pipeline identified VCP as a drug target for frontotemporal dementia, how would CAR-T approaches compare to small molecule inhibitors for this target?"

**Talking points:**

- "VCP is the demo target from Stage 2 — identified by the genomics pipeline as a frontotemporal dementia drug target."
- "The CAR-T agent evaluates whether cell therapy approaches could work alongside small molecule inhibitors."
- "It pulls from literature on CAR-T for neurodegeneration, constructs targeting intracellular proteins, and manufacturing considerations."
- "This shows the power of cross-functional intelligence — one agent informing another's decision space."

---

### Cross-Functional Comparative Query (3 minutes)

```bash
curl -s -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Compare manufacturing parameters and clinical outcomes for CD19 vs BCMA CAR-T products"
  }' | python3 -m json.tool
```

**Talking points:**

- "Comparative analysis across manufacturing AND clinical collections simultaneously."
- "CD19 products: 4 approved (Kymriah, Yescarta, Tecartus, Breyanzi) — compare transduction efficiency, expansion time, release criteria."
- "BCMA products: 2 approved (Abecma, Carvykti) — different manufacturing challenges, different patient populations."
- "This cross-functional view doesn't exist in any single publication — it's synthesized from 11 collections."

---

### Stage 3 Integration (3 minutes)

**Show:** Drug Discovery UI at http://localhost:8505

**Talking points:**

- "When the CAR-T agent identifies resistance mechanisms — like BCMA downregulation in myeloma — those targets feed into Stage 3."
- "BioNeMo generates small molecule candidates that could be combined with CAR-T therapy."
- "Combination strategies: CAR-T for bulk disease + small molecule for resistant clones."
- "The platform enables multi-modal therapeutic design — not just cell therapy OR drug discovery, but both together."

---

### Imaging Agent Connection (3 minutes)

**Talking points:**

- "The Imaging Intelligence Agent runs on the same platform. When it detects a high-risk lung nodule, cross-modal triggers query the genomic evidence."
- "If the genomic analysis identifies a target relevant to both solid tumor CAR-T and small molecule inhibition, both agents contribute intelligence."
- "Example: Imaging finds lung mass → genomics identifies EGFR → CAR-T agent evaluates EGFR-targeted cell therapy → Drug Discovery generates EGFR inhibitors."
- "Multi-agent collaboration, all sharing the same 3.56 million variant truth."

---

### The Complete Loop (3 minutes)

**Talking points:**

Walk through the full precision medicine loop:

```
Patient DNA (Illumina Sequencing)
    ↓
Stage 1: GPU Genomics (Parabricks)
    ↓ 11.7M variants called
Stage 2: RAG Target Identification (Claude + Milvus)
    ↓ VCP identified as drug target
    ↓ genomic_evidence: 3.56M vectors (shared)
    ├─→ CAR-T Intelligence Agent
    │   └─ Evaluates cell therapy approaches for target
    │   └─ Cross-references 6,266 therapy-specific vectors
    │   └─ Identifies optimal construct design + manufacturing
    ├─→ Imaging Intelligence Agent
    │   └─ Cross-modal triggers from imaging findings
    │   └─ Connects phenotype to genotype
    └─→ Stage 3: Drug Discovery (BioNeMo)
        └─ 100 novel drug candidates generated
        └─ Docking + drug-likeness scoring

Combined Clinical Output:
    → FHIR R4 Reports (imaging)
    → PDF Reports (CAR-T + drug discovery)
    → JSON data for dashboards
```

- "Patient DNA to therapeutic strategy. Genomics, cell therapy intelligence, imaging AI, and drug discovery — all on one platform."
- "Every agent sees the same genomic truth. Every answer is grounded in evidence."
- "This is the future of precision medicine — multi-modal, multi-agent, evidence-grounded."

---

### Closing Route B (2 minutes)

**Talking points:**

- "The CAR-T Intelligence Agent isn't standalone — it's part of a precision medicine ecosystem."
- "Shared infrastructure means shared intelligence. 3.56 million genomic vectors inform every agent."
- "From a $3,999 DGX Spark to enterprise-scale DGX SuperPOD — the same code, the same agents, the same evidence."
- "All Apache 2.0. Build on it, extend it, deploy it."

**Scaling story:**

| Phase | Hardware | Scale |
|---|---|---|
| Phase 1 | DGX Spark ($3,999) | Proof build — what you just saw |
| Phase 2 | DGX B200 | Department — multi-patient analysis |
| Phase 3 | DGX SuperPOD | Enterprise — federated multi-site intelligence |

---

## Troubleshooting

### Milvus Connection Issues

```bash
# Check Milvus status
curl -s http://localhost:19530/v1/vector/collections

# Check dependencies
docker compose logs milvus-etcd
docker compose logs milvus-minio

# If collections empty, re-run setup
python3 scripts/setup_collections.py --seed-constructs
python3 scripts/ingest_pubmed.py --max-results 5000
python3 scripts/ingest_clinical_trials.py --max-results 1500
```

### Claude API Timeout

Queries take ~24 seconds due to LLM generation. If timeouts occur:

```bash
# Verify API key
echo $ANTHROPIC_API_KEY | head -c 10

# Test Claude directly
curl -s https://api.anthropic.com/v1/messages \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "anthropic-version: 2023-06-01" \
  -H "content-type: application/json" \
  -d '{"model": "claude-sonnet-4-20250514", "max_tokens": 10, "messages": [{"role": "user", "content": "Hi"}]}'
```

### Comparative Mode Not Triggering

Comparative mode requires keywords: "compare", "vs", "versus", or "comparing". Ensure one of these words appears in the query.

### Empty Search Results

If searches return no results, verify collection data:

```bash
curl -s http://localhost:8522/health | python3 -m json.tool
# Check that cart_literature shows 5047, cart_trials shows 973
```

If counts are zero, re-run ingestion scripts (see Pre-Demo Setup).

### PDF Export Issues

PDF generation requires ReportLab. If PDF export fails:

```bash
pip install reportlab
```

---

## Quick Reference

| Action | Command / URL |
|---|---|
| Health check | `curl http://localhost:8522/health` |
| Knowledge graph | `curl http://localhost:8522/knowledge/stats` |
| List collections | `curl http://localhost:8522/collections` |
| RAG query | `curl -X POST http://localhost:8522/api/ask` |
| Full RAG query | `curl -X POST http://localhost:8522/query` |
| Evidence search | `curl -X POST http://localhost:8522/search` |
| Find related | `curl -X POST http://localhost:8522/find-related` |
| Markdown report | `curl -X POST http://localhost:8522/api/report/markdown` |
| JSON report | `curl -X POST http://localhost:8522/api/report/json` |
| PDF report | `curl -X POST http://localhost:8522/api/report/pdf` |
| CAR-T Agent UI | http://localhost:8521 |
| CAR-T Agent API docs | http://localhost:8522/docs |
| RAG Chat (Stage 2) | http://localhost:8501 |
| Drug Discovery UI | http://localhost:8505 |
| Landing Page | http://localhost:8080 |
| Grafana Monitoring | http://localhost:3000 |
| Prometheus Metrics | `curl http://localhost:8522/metrics` |

---

*HCLS AI Factory — Apache 2.0 | February 2026*
