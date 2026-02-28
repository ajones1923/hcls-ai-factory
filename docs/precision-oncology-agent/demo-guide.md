---
search:
  boost: 2
tags:
  - Demo
  - Walkthrough
  - Precision Oncology
  - Genomics
  - MTB
  - Clinical Decision Support
---

# Precision Oncology Agent — Demo Guide

> **Step-by-step walkthrough for demonstrating the Precision Oncology Agent on DGX Spark.**
>
> License: Apache 2.0 | Date: February 2026

---

## Demo Overview

| Parameter | Value |
|---|---|
| **Route A Duration** | 20 minutes (standalone agent) |
| **Route B Duration** | 30 minutes (cross-platform integration) |
| **Hardware** | NVIDIA DGX Spark (GB10, 128 GB unified) |
| **Knowledge Base** | ~1,490 vectors across 10 owned collections + 3.5M genomic vectors |
| **Knowledge Graph** | ~40 actionable targets, ~30 therapies, ~20 resistance mechanisms |
| **LLM** | Claude Sonnet 4.6 (Anthropic) |
| **Export Formats** | Markdown, JSON, PDF, FHIR R4 |

### What the Audience Will See

1. A health dashboard showing 11 vector collections spanning the precision oncology clinical workflow
2. Cross-collection RAG queries pulling evidence from variants, literature, therapies, guidelines, trials, biomarkers, resistance, pathways, outcomes, and cases -- simultaneously
3. Automatic comparative analysis triggered by natural language -- "Compare X vs Y" produces structured tables with efficacy, safety, and guideline data
4. MTB packet generation from raw VCF data -- variant annotation, therapy ranking, trial matching, and open questions in under 30 seconds
5. Evidence-based therapy ranking with AMP/ASCO/CAP evidence tiers (A-D) and resistance awareness
6. Hybrid clinical trial matching combining deterministic filters with semantic search and composite scoring
7. Clickable PubMed and ClinicalTrials.gov citations grounding every answer
8. Professional reports exported as Markdown, JSON, NVIDIA-themed PDF, and FHIR R4 DiagnosticReport Bundles
9. *(Route B)* Cross-pipeline data sharing -- oncology intelligence layered on 3.5 million genomic variant vectors from Stage 1
10. *(Route B)* The complete loop: Patient DNA -> Variant Analysis -> Precision Oncology Intelligence -> Drug Candidates

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
cd ai_agent_adds/precision_oncology_agent
docker compose up -d

# This starts:
# Milvus (etcd + MinIO + standalone)
# Streamlit MTB Workbench (port 8526)
# FastAPI server (port 8527)
```

### Step 4: Verify All Services Healthy

```bash
curl -s http://localhost:8527/health | python3 -m json.tool
```

Expected response:

```json
{
  "status": "healthy",
  "collections": {
    "onco_variants": 300,
    "onco_literature": 500,
    "onco_therapies": 120,
    "onco_guidelines": 100,
    "onco_trials": 200,
    "onco_biomarkers": 80,
    "onco_resistance": 80,
    "onco_pathways": 50,
    "onco_outcomes": 50,
    "onco_cases": 10,
    "genomic_evidence": 3561170
  },
  "total_vectors": 3562660,
  "version": "0.1.0",
  "services": {
    "milvus": true,
    "embedder": true,
    "rag_engine": true,
    "intelligence_agent": true,
    "case_manager": true,
    "trial_matcher": true,
    "therapy_ranker": true
  }
}
```

All 11 collections should have non-zero counts. All 7 services should be `true`.

### Step 5: Verify Knowledge Graph

```bash
curl -s http://localhost:8527/knowledge/stats | python3 -m json.tool
```

Expected: ~40 actionable targets, ~120 therapies, ~80 resistance mechanisms, ~50 pathways, ~80 biomarkers.

---

## Route A: Standalone Agent Demo (20 minutes)

### Step 1: Opening (1 minute)

**Talking points:**

- "This is the Precision Oncology Agent -- it generates Molecular Tumor Board packets from raw genomic data."
- "11 Milvus collections covering variants, literature, therapies, guidelines, trials, biomarkers, resistance mechanisms, signaling pathways, treatment outcomes, and patient cases."
- "Every query searches all data sources simultaneously. Claude synthesizes cross-functional clinical insights with citations."

**Show:** Health endpoint -- highlight 11 collections, ~1,490 owned vectors, 3.5M genomic vectors, 7 services active.

---

### Step 2: RAG Knowledge Query (3 minutes)

```bash
curl -s -X POST http://localhost:8527/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What therapies target BRAF V600E in melanoma?",
    "cancer_type": "melanoma",
    "gene": "BRAF"
  }' | python3 -m json.tool
```

Expected response highlights:

```json
{
  "answer": "BRAF V600E melanoma has several FDA-approved targeted therapy options...",
  "sources": [
    {"collection": "onco_therapies", "score": 0.91, "text": "Vemurafenib (Zelboraf) is a selective BRAF V600E inhibitor..."},
    {"collection": "onco_variants", "score": 0.88, "text": "BRAF V600E is the most common BRAF mutation in melanoma..."},
    {"collection": "onco_guidelines", "score": 0.85, "text": "NCCN Melanoma v2.2025: BRAF-targeted therapy recommended..."},
    {"collection": "onco_resistance", "score": 0.82, "text": "BRAF inhibitor resistance via NRAS activation or MEK1/2 mutations..."}
  ],
  "follow_up_questions": [
    "What resistance mechanisms emerge after BRAF inhibitor therapy?",
    "Should BRAF+MEK combination be used over BRAF monotherapy?",
    "What is the role of immunotherapy in BRAF-mutant melanoma?"
  ],
  "confidence": 0.89,
  "processing_time_ms": 24100
}
```

**Talking points:**

- "One question, four collections hit simultaneously -- therapies, variants, guidelines, and resistance."
- "The answer identifies specific drugs: vemurafenib, dabrafenib+trametinib combination, and encorafenib+binimetinib."
- "Every claim is backed by a citation with cosine similarity scores. The top hit scored 0.91."
- "Claude suggests three follow-up questions, each targeting a different clinical dimension."

---

### Step 3: Comparative Analysis (3 minutes)

```bash
curl -s -X POST http://localhost:8527/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Compare osimertinib vs erlotinib for EGFR-mutant NSCLC"
  }' | python3 -m json.tool
```

**Talking points:**

- "Watch what happens when I say 'compare'. The engine automatically detects this is a comparative query."
- "It parses two entities -- osimertinib and erlotinib -- and resolves each against the knowledge graph."
- "Dual retrieval runs: entity A and entity B are searched separately, then results are merged."
- "Claude produces an 8-section comparison: MoA differences, efficacy data, safety profile, biomarker considerations, resistance mechanisms, guideline recommendations, trial evidence, and summary."

**Expected output structure:**

```
## Comparison: Osimertinib vs Erlotinib

| Dimension | Osimertinib (Tagrisso) | Erlotinib (Tarceva) |
|---|---|---|
| Generation | 3rd-generation EGFR TKI | 1st-generation EGFR TKI |
| Targets | EGFR exon19del/L858R + T790M | EGFR exon19del/L858R |
| CNS penetration | High | Limited |
| Median PFS (first-line) | 18.9 months (FLAURA) | 10.2 months (FLAURA) |
| T790M resistance | Active against | No activity |
| Key trial | FLAURA, AURA3 | EURTAC, OPTIMAL |

### NCCN Recommendation
Osimertinib is preferred first-line for EGFR-mutant NSCLC (Category 1).
```

**Comparison types you can demo:**

| Query | Entities Resolved |
|---|---|
| "Compare osimertinib vs erlotinib" | Drug vs drug |
| "BRAF+MEK inhibition vs immunotherapy for melanoma" | Modality comparison |
| "Pembrolizumab versus nivolumab" | Product vs product |
| "Compare PARP inhibitors for BRCA-mutant ovarian cancer" | Drug class comparison |

---

### Step 4: Create a Patient Case (3 minutes)

```bash
curl -s -X POST http://localhost:8527/api/cases \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "ONCO-DEMO-001",
    "cancer_type": "NSCLC",
    "stage": "IV",
    "variants": [
      {"gene": "EGFR", "variant": "L858R", "variant_type": "SNV"},
      {"gene": "TP53", "variant": "R273H", "variant_type": "SNV"},
      {"gene": "KRAS", "variant": "G12C", "variant_type": "SNV"}
    ],
    "biomarkers": {
      "MSI": "MSS",
      "TMB": 14.2,
      "PD-L1_TPS": 80
    },
    "prior_therapies": ["carboplatin", "pemetrexed"]
  }' | python3 -m json.tool
```

Expected response:

```json
{
  "case_id": "a1b2c3d4-e5f6-7890-abcd-ef1234567890",
  "patient_id": "ONCO-DEMO-001",
  "cancer_type": "NSCLC",
  "stage": "IV",
  "variant_count": 3,
  "biomarkers": {"MSI": "MSS", "TMB": 14.2, "PD-L1_TPS": 80},
  "prior_therapies": ["carboplatin", "pemetrexed"],
  "created_at": "2026-02-28T12:00:00Z"
}
```

**Talking points:**

- "We just created a patient case with 3 variants: EGFR L858R, TP53 R273H, and KRAS G12C."
- "Biomarkers include TMB at 14.2 (TMB-high), PD-L1 TPS at 80%, and MSI-stable."
- "Prior therapies are carboplatin and pemetrexed -- the agent will check for resistance patterns."
- "The case is stored in the onco_cases collection for MTB packet generation."

---

### Step 5: Generate MTB Packet (4 minutes)

```bash
# Use the case_id from the previous step
curl -s -X POST http://localhost:8527/api/cases/{CASE_ID}/mtb \
  -H "Content-Type: application/json" \
  -d '{
    "include_trials": true,
    "include_therapies": true,
    "include_resistance": true,
    "top_k": 10
  }' | python3 -m json.tool
```

**Talking points:**

- "This is the core product -- a complete Molecular Tumor Board packet generated in under 30 seconds."
- "The variant table shows 3 variants with actionability classification:"
  - "EGFR L858R: Level A -- FDA-approved companion diagnostic exists. Drugs: osimertinib, erlotinib, gefitinib."
  - "KRAS G12C: Level A -- sotorasib (Lumakras) and adagrasib (Krazati) are FDA-approved."
  - "TP53 R273H: VUS -- no direct actionable therapy, but flags tumor suppressor loss."
- "The therapy ranking places osimertinib at #1 (Level A evidence, NCCN Category 1) and sotorasib at #2."
- "Resistance check: carboplatin and pemetrexed in prior therapies -- the agent flags no cross-resistance with osimertinib."
- "Trial matches: 10 recruiting trials matched by cancer type, biomarker criteria, and stage."
- "Open questions flag the TP53 VUS and suggest functional assay consideration."

---

### Step 6: Trial Matching (2 minutes)

```bash
curl -s -X POST http://localhost:8527/api/trials/match \
  -H "Content-Type: application/json" \
  -d '{
    "cancer_type": "NSCLC",
    "biomarkers": {
      "EGFR": "L858R",
      "PD-L1_TPS": 80,
      "TMB": 14.2
    },
    "stage": "IV",
    "top_k": 5
  }' | python3 -m json.tool
```

**Talking points:**

- "Hybrid trial matching: deterministic filter on cancer type + recruiting status, plus semantic search on the full patient profile."
- "Composite scoring: 40% biomarker match, 25% semantic similarity, 20% phase weight, 15% status weight."
- "Each match includes: trial ID, phase, match score, matched criteria, unmatched criteria, and explanation."
- "The top match is a Phase 3 EGFR-mutant NSCLC trial with 0.87 composite score -- all biomarkers matched."

---

### Step 7: Therapy Ranking (2 minutes)

```bash
curl -s -X POST http://localhost:8527/api/therapies/rank \
  -H "Content-Type: application/json" \
  -d '{
    "cancer_type": "NSCLC",
    "variants": [
      {"gene": "EGFR", "variant": "L858R"},
      {"gene": "KRAS", "variant": "G12C"}
    ],
    "biomarkers": {
      "TMB": 14.2,
      "PD-L1_TPS": 80
    },
    "prior_therapies": ["carboplatin", "pemetrexed"]
  }' | python3 -m json.tool
```

**Talking points:**

- "Six-step therapy ranking: variant-driven identification, biomarker-driven identification, evidence-level sort, resistance check, contraindication check, and evidence retrieval."
- "Variant-driven: EGFR L858R -> osimertinib (Level A); KRAS G12C -> sotorasib (Level A)."
- "Biomarker-driven: TMB-H (14.2) -> pembrolizumab (Level A); PD-L1 TPS 80% -> pembrolizumab first-line (Level A)."
- "No resistance flags raised against prior carboplatin/pemetrexed."
- "Each therapy has supporting evidence from onco_therapies and onco_literature with citations."

---

### Step 8: Report Export (2 minutes)

```bash
# Markdown report
curl -s -X POST http://localhost:8527/api/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the optimal treatment sequence for EGFR-mutant NSCLC?",
    "cancer_type": "NSCLC",
    "format": "markdown"
  }' | head -50

# JSON report
curl -s -X POST http://localhost:8527/api/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the optimal treatment sequence for EGFR-mutant NSCLC?",
    "format": "json"
  }' | python3 -m json.tool | head -30

# PDF report
curl -s -X POST http://localhost:8527/api/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the optimal treatment sequence for EGFR-mutant NSCLC?",
    "format": "pdf"
  }' --output nsclc_egfr_report.pdf

# FHIR R4 export (for existing case)
curl -s http://localhost:8527/api/reports/{CASE_ID}/fhir | python3 -m json.tool | head -30
```

**Talking points:**

- "Four export formats: Markdown for sharing, JSON for programmatic consumption, PDF for clinical documentation, FHIR R4 for interoperability."
- "The PDF is NVIDIA-themed with professional styling via ReportLab."
- "The FHIR R4 export produces a DiagnosticReport Bundle coded with SNOMED CT and LOINC -- ready for EHR integration."
- "Every report includes: query metadata, analysis, evidence sources, knowledge graph context, and citation links."

---

### Step 9: Closing Route A (1 minute)

**Talking points:**

- "11 collections, ~1,490 curated vectors, ~40 actionable targets, 6-step therapy ranking, 4-step trial matching."
- "This is a complete precision oncology clinical decision support system -- from VCF to MTB packet."
- "Every answer is grounded in published evidence with clickable citations. No hallucination."
- "AMP/ASCO/CAP evidence tiers ensure clinicians know the strength of each recommendation."

---

## Route B: Cross-Platform Integration Demo (30 minutes)

> **Prerequisite:** Complete Route A first, or start from a fresh session with all services running.
>
> This route demonstrates how the Precision Oncology Agent connects to the full HCLS AI Factory -- bridging genomic variant analysis, clinical decision support, imaging intelligence, and AI-driven drug discovery.

### Step 1: Platform Overview (2 minutes)

**Show:** HCLS AI Factory landing page at http://localhost:8080

**Talking points:**

- "The HCLS AI Factory is a 3-stage precision medicine platform: GPU-accelerated genomics, RAG-grounded target identification, and AI-driven drug discovery."
- "The Precision Oncology Agent extends this platform with clinical decision support for Molecular Tumor Boards."
- "All agents share the same Milvus vector database, BGE embeddings, and Claude LLM -- creating a unified intelligence layer."

**Architecture overview:**

```
Stage 1: GPU Genomics (Parabricks)
    | 11.7M variants -> 3.5M quality-filtered
Stage 2: RAG Target Identification (Claude + Milvus)
    | genomic_evidence: 3.5M vectors (shared)
Precision Oncology Agent (11 collections)
    | VCF -> Variant Annotation -> Therapy Ranking -> MTB Packet
Stage 3: Drug Discovery (BioNeMo)
    | 100 ranked drug candidates
```

---

### Step 2: Shared Genomic Data Layer (3 minutes)

```bash
# Show the Oncology agent's view of all collections
curl -s http://localhost:8527/collections | python3 -m json.tool
```

**Talking points:**

- "The Oncology agent has 10 specialized collections -- ~1,490 curated vectors covering variants, literature, therapies, guidelines, trials, biomarkers, resistance, pathways, outcomes, and cases."
- "It also reads from `genomic_evidence` -- 3,561,170 vectors from the genomics pipeline. These are the same variant annotations Stage 2 uses."
- "This is the key integration point. The same patient variants that Stage 2 analyzes for drug targets are available to inform oncology treatment decisions."

---

### Step 3: Genomic Evidence Integration (3 minutes)

```bash
curl -s -X POST http://localhost:8527/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What genomic variants in this patient affect targeted therapy selection and resistance prediction?",
    "cancer_type": "NSCLC",
    "gene": "EGFR"
  }' | python3 -m json.tool
```

**Talking points:**

- "This query bridges genomic variants and oncology treatment intelligence. Claude searches both the therapy collections and the genomic evidence."
- "It identifies EGFR mutations relevant to TKI selection -- L858R, exon 19 deletions, T790M resistance."
- "The genomic_evidence collection provides variant-level annotations from ClinVar and AlphaMissense."
- "This is precision medicine at work -- the patient's genome informs the treatment recommendation."

---

### Step 4: End-to-End MTB Workflow (5 minutes)

**Show:** Streamlit MTB Workbench at http://localhost:8526

**Walk through the 5-tab workflow:**

1. **Case Management Tab:** Create a case with cancer type NSCLC, stage IV, EGFR L858R variant, TMB 14.2, PD-L1 TPS 80%.
2. **Evidence Explorer Tab:** Ask "What is the optimal first-line treatment for EGFR L858R NSCLC?" -- show streaming answer with citations.
3. **Trial Matching Tab:** Match the case to recruiting trials -- show top 5 matches with composite scores.
4. **Therapy Ranking Tab:** View ranked therapies with resistance flags -- osimertinib at #1.
5. **Outcomes Dashboard Tab:** Export the MTB packet as PDF -- download and open to show the NVIDIA-themed report.

**Talking points:**

- "Five tabs, one integrated workflow. From case creation to MTB packet in under 2 minutes."
- "The evidence explorer provides interactive clinical Q&A with conversation memory."
- "Trial matching runs in real-time against 200+ oncology trials."
- "The therapy ranker integrates variant-level, biomarker-level, and resistance-aware analysis."

---

### Step 5: Cross-Agent Collaboration (3 minutes)

```bash
curl -s -X POST http://localhost:8527/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "If imaging detects a new lung lesion in a patient with EGFR-mutant NSCLC on osimertinib, what resistance mechanisms should be evaluated and what treatment options exist?"
  }' | python3 -m json.tool
```

**Talking points:**

- "This query spans multiple domains: imaging findings, resistance biology, and treatment options."
- "The cross-modal trigger fires when Level A/B actionable variants are detected -- it queries both genomic evidence and imaging collections."
- "The agent identifies EGFR C797S and MET amplification as the most common osimertinib resistance mechanisms."
- "It suggests: amivantamab-vmjw for EGFR/MET bispecific targeting, or combination with MET inhibitors like savolitinib."
- "This is multi-agent intelligence -- imaging findings trigger genomic re-analysis, which informs treatment adaptation."

---

### Step 6: Drug Discovery Pipeline Connection (3 minutes)

**Show:** Drug Discovery UI at http://localhost:8505

**Talking points:**

- "When the Oncology agent identifies actionable targets, those targets feed into Stage 3."
- "BioNeMo generates small molecule candidates via MolMIM and evaluates binding via DiffDock."
- "Example: EGFR T790M resistance identified -> drug discovery pipeline generates novel T790M-selective molecules."
- "Combination strategies: osimertinib for primary EGFR mutation + novel molecule for resistance clone."
- "The platform enables rational drug design informed by clinical resistance data."

---

### Step 7: The Complete Loop (3 minutes)

**Talking points:**

Walk through the full precision medicine loop:

```
Patient DNA (Illumina Sequencing)
    |
Stage 1: GPU Genomics (Parabricks)
    | 11.7M variants called
Stage 2: RAG Target Identification (Claude + Milvus)
    | genomic_evidence: 3.5M vectors (shared)
    |
    +---> Precision Oncology Agent
    |     +-- VCF parsed, variants annotated
    |     +-- EGFR L858R: Level A, osimertinib recommended
    |     +-- TMB-H + PD-L1: pembrolizumab eligible
    |     +-- 8 recruiting trials matched
    |     +-- MTB packet exported (PDF + FHIR R4)
    |
    +---> Imaging Intelligence Agent
    |     +-- Cross-modal triggers from imaging findings
    |     +-- Connects phenotype to genotype
    |
    +---> Stage 3: Drug Discovery (BioNeMo)
          +-- 100 novel drug candidates generated
          +-- Docking + drug-likeness scoring

Combined Clinical Output:
    -> MTB Packet (PDF, FHIR R4) from Precision Oncology
    -> Imaging Reports (FHIR R4) from Imaging Intelligence
    -> Drug Candidate Report from Drug Discovery
```

- "Patient DNA to therapeutic strategy. Genomics, oncology intelligence, imaging AI, and drug discovery -- all on one platform."
- "Every agent sees the same genomic truth. Every answer is grounded in evidence."
- "Interoperability through FHIR R4 -- SNOMED CT and LOINC coded for EHR integration."

---

### Closing Route B (2 minutes)

**Talking points:**

- "The Precision Oncology Agent is not standalone -- it is part of a precision medicine ecosystem."
- "Shared infrastructure means shared intelligence. 3.5 million genomic vectors inform every agent."
- "AMP/ASCO/CAP evidence tiers provide the clinical rigor oncologists require."
- "From a $3,999 DGX Spark to enterprise-scale DGX SuperPOD -- the same code, the same agents, the same evidence."
- "All Apache 2.0. Build on it, extend it, deploy it."

**Scaling story:**

| Phase | Hardware | Scale |
|---|---|---|
| Phase 1 | DGX Spark ($3,999) | Proof build -- what you just saw |
| Phase 2 | DGX B200 | Department -- multi-patient MTB workflow |
| Phase 3 | DGX SuperPOD | Enterprise -- federated multi-site tumor boards |

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
python3 scripts/setup_collections.py --seed
python3 scripts/ingest_pubmed.py --max-results 5000
python3 scripts/ingest_clinical_trials.py --max-results 1500
python3 scripts/ingest_civic.py
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

Comparative mode requires keywords: "compare", "vs", "versus", "difference between", or "head-to-head". Ensure one of these words appears in the query.

### Empty Search Results

If searches return no results, verify collection data:

```bash
curl -s http://localhost:8527/health | python3 -m json.tool
# Check that onco_variants shows ~300, onco_literature shows ~500
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
| Health check | `curl http://localhost:8527/health` |
| Knowledge graph | `curl http://localhost:8527/knowledge/stats` |
| List collections | `curl http://localhost:8527/collections` |
| RAG query | `curl -X POST http://localhost:8527/api/ask` |
| Full RAG query | `curl -X POST http://localhost:8527/query` |
| Evidence search | `curl -X POST http://localhost:8527/search` |
| Find related | `curl -X POST http://localhost:8527/find-related` |
| Create case | `curl -X POST http://localhost:8527/api/cases` |
| Generate MTB packet | `curl -X POST http://localhost:8527/api/cases/{id}/mtb` |
| Match trials | `curl -X POST http://localhost:8527/api/trials/match` |
| Rank therapies | `curl -X POST http://localhost:8527/api/therapies/rank` |
| Generate report | `curl -X POST http://localhost:8527/api/reports/generate` |
| Export case (markdown) | `curl http://localhost:8527/api/reports/{id}/markdown` |
| Export case (json) | `curl http://localhost:8527/api/reports/{id}/json` |
| Export case (pdf) | `curl http://localhost:8527/api/reports/{id}/pdf` |
| Export case (fhir) | `curl http://localhost:8527/api/reports/{id}/fhir` |
| Oncology Agent UI | http://localhost:8526 |
| Oncology Agent API docs | http://localhost:8527/docs |
| RAG Chat (Stage 2) | http://localhost:8501 |
| Drug Discovery UI | http://localhost:8505 |
| Landing Page | http://localhost:8080 |
| Grafana Monitoring | http://localhost:3000 |
| Prometheus Metrics | `curl http://localhost:8527/metrics` |

---

*HCLS AI Factory -- Apache 2.0 | February 2026*
