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
| **Route B Duration** | 25 minutes (cross-platform integration) |
| **Hardware** | NVIDIA DGX Spark (GB10, 128 GB unified) |
| **Knowledge Base** | 6,266+ vectors across 10 owned collections + 3.56M genomic vectors |
| **Knowledge Graph** | 25 target antigens, 6 FDA products, 8 toxicity profiles |
| **LLM** | Claude Sonnet 4.6 (Anthropic) |
| **Export Formats** | Markdown, JSON, PDF |
| **Primary UI** | Streamlit at http://localhost:8521 |
| **API** | FastAPI at http://localhost:8522 |
| **Tests** | 241 tests |

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
10. *(Route B)* The complete loop: Patient DNA -> Variant Analysis -> CAR-T Intelligence -> Drug Candidates

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

### Step 6: Open Streamlit UI

Open http://localhost:8521 in a browser. Confirm the Chat tab loads with the sidebar showing all 10 collection checkboxes checked by default.

---

## Route A: Standalone Agent Demo (15 minutes)

> **All interaction through the Streamlit UI at http://localhost:8521. Zero curl commands, zero terminal windows.**

| Step | Tab / UI | Action | Time |
|------|----------|--------|------|
| 1 — Opening | Landing page `:8080` | Show health grid, CAR-T service green | 1 min |
| 2 — Evidence Query | Chat tab `:8521` | Type query about CD19 CAR-T failure | 3 min |
| 3 — Sidebar Filters | Chat tab sidebar | Show 10 collection checkboxes, select target, demo query buttons | 1 min |
| 4 — Comparative | Chat tab `:8521` | Type "Compare 4-1BB vs CD28 costimulatory domains for DLBCL" | 3 min |
| 5 — Knowledge Graph | Knowledge Graph tab | Select entity type, explore graph visualization | 2 min |
| 6 — Deep Research | Chat tab (toggle) | Enable "Deep Research Mode" in sidebar, type deep research query | 2 min |
| 7 — Export | Chat tab | Click "Download PDF" button after query results appear | 1 min |
| 8 — Closing | -- | Talking points: 10 collections, 6,266+ vectors, grounded citations | 2 min |

---

### Step 1: Opening (1 minute)

**Show:** Landing page at http://localhost:8080 -- confirm the CAR-T Intelligence service shows green in the health grid.

**Talking points:**

- "This is the CAR-T Intelligence Agent — it breaks down data silos across the entire CAR-T cell therapy development lifecycle."
- "11 Milvus collections covering literature, clinical trials, FDA-approved constructs, assays, manufacturing, safety, biomarkers, regulatory, sequences, real-world evidence, and genomic variants."
- "Every query searches all data sources simultaneously and Claude synthesizes cross-functional insights with citations."

**Show:** CAR-T Agent UI at http://localhost:8521. Point out:

- The header: "CAR-T Intelligence Agent" in NVIDIA green with the subtitle "Cross-Functional Intelligence Across the CAR-T Development Lifecycle"
- Three tabs across the top: **Chat**, **Knowledge Graph**, **Image Analysis**
- The sidebar with Configuration section, collection checkboxes with live vector counts, and Demo Queries section

---

### Step 2: RAG Knowledge Query (3 minutes)

**Show:** Chat tab at http://localhost:8521

**Select:** In the sidebar, verify all 10 collection checkboxes are checked. Each checkbox shows the collection name and live record count (e.g., "Literature (5,047)", "Clinical Trials (973)").

**Select:** Target Antigen Filter = `CD19`

**Type this query:**

> Why do CD19 CAR-T therapies fail in relapsed B-ALL?

**Expected result:** A status indicator appears: "Searching across CAR-T data sources..." with real-time progress. Claude returns a multi-paragraph answer streamed token-by-token. Look for:

- An expandable **"Evidence Sources"** panel showing hit count and search time (e.g., "Evidence Sources (24 results, 14ms)")
- Inside the evidence panel: color-coded evidence cards from multiple collections:
  - **cart_literature** evidence card (blue badge, score ~0.89) — "Antigen loss observed in 28% of relapses..."
  - **cart_assays** evidence card (yellow badge, score ~0.85) — "CD19-negative relapses via trogocytosis..."
  - **cart_trials** evidence card (purple badge, score ~0.82) — "ELIANA trial: 12-month EFS 73%..."
  - **cart_constructs** evidence card (green badge, score ~0.80) — "Kymriah tisagenlecleucel 4-1BB..."
- Clickable PubMed and ClinicalTrials.gov links on relevant evidence cards
- Three follow-up question suggestions below the answer
- Three export buttons below the response: **"Download Markdown"**, **"Download JSON"**, **"Download PDF"**

**Talking points:**

- "One question, four collections hit simultaneously — literature, assays, trials, and constructs."
- "The answer identifies specific failure mechanisms: antigen loss in 28% of relapses, lineage switch in 10% with KMT2A rearrangements, trogocytosis, and T-cell exhaustion."
- "Every claim is backed by a citation with cosine similarity scores. The top hit scored 0.89 — very high relevance."
- "Evidence cards are color-coded by collection — blue for literature, purple for trials, green for constructs, yellow for assays."
- "Claude suggests three follow-up questions, each targeting a different aspect of the problem."

---

### Step 3: Sidebar Filters and Demo Queries (1 minute)

**Show:** The sidebar at http://localhost:8521

- Point out the **10 collection checkboxes** — each one maps to a Milvus collection with a live record count displayed.
- Point out the **Target Antigen Filter** selectbox (All Targets, CD19, BCMA, CD22, CD20, CD30, CD33, CD38, CD123, GD2, HER2, GPC3, EGFR, Mesothelin, PSMA, ROR1).
- Point out the **Development Stage** selectbox (All Stages, Target Identification, CAR Design, Vector Engineering, Testing, Clinical).
- Point out the **Date Range** section — two year inputs (From Year / To Year) and an "Apply date filter" checkbox.
- Point out the **Total** line at the bottom of the collections section showing total vector count across selected collections.

**Click:** One of the 13 pre-built **Demo Query** buttons in the sidebar (e.g., "BCMA CAR-T resistance mechanisms in myeloma").

**Expected result:** The demo query auto-populates the chat input and executes. A fully cited RAG answer appears with evidence cards and export buttons.

**Talking points:**

- "The sidebar gives full control over which evidence sources are queried. Uncheck collections to narrow the scope."
- "Pre-built demo queries let clinicians start with validated questions — no prompt engineering required."
- "13 demo queries cover the full development lifecycle — from construct binding affinity to real-world outcomes to regulatory pathways."

---

### Step 4: Comparative Analysis (3 minutes)

**Show:** Chat tab at http://localhost:8521

**Type this query:**

> Compare 4-1BB vs CD28 costimulatory domains for DLBCL

**Expected result:** The engine detects "Compare" and triggers dual-entity retrieval. The status indicator updates: "Comparative analysis: **4-1BB** vs **CD28**" with hit count and search time. Claude produces a structured comparison:

- An expandable **"Comparative Evidence"** panel with two sections:
  - Entity A header (blue): **4-1BB** — evidence cards for 4-1BB
  - A "-- VS --" divider in NVIDIA green
  - Entity B header (purple): **CD28** — evidence cards for CD28
- Claude's synthesized comparison table:

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

**Talking points:**

- "Watch what happens when I say 'compare'. The engine automatically detects this is a comparative query."
- "It parses two entities — 4-1BB and CD28 — and resolves each against the knowledge graph."
- "Dual retrieval runs: entity A and entity B are searched separately, then results are merged into a structured prompt."
- "The evidence panel shows both sides — 4-1BB hits in blue, CD28 hits in purple, separated by a VS divider."
- "Claude produces a comparison table with advantages, limitations, and clinical context."
- "Five of six FDA-approved CAR-T products use 4-1BB. Only Yescarta uses CD28. That's a clear trend."

**Comparison types you can demo:**

| Query | Entities Resolved |
|---|---|
| "Compare CD19 vs BCMA" | Target antigen vs target antigen |
| "Kymriah versus Carvykti" | Product -> CD19 vs BCMA |
| "Compare CRS and ICANS" | Toxicity vs toxicity |
| "4-1BB vs CD28 for persistence" | Costimulatory domain comparison |
| "Lentiviral vs transposon vectors" | Manufacturing comparison |

---

### Step 5: Knowledge Graph Exploration (2 minutes)

**Click:** The **"Knowledge Graph"** tab at the top of the page at http://localhost:8521.

**Select:** Entity Type = `Target Antigens` from the selectbox on the left.

**Expected result:** An interactive pyvis graph visualization renders showing:

- Target antigen nodes in NVIDIA green (CD19, BCMA, CD22, etc.)
- Disease nodes in blue (B-ALL, DLBCL, FL, MCL, Multiple Myeloma)
- FDA-approved product nodes in purple (Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti)
- Resistance mechanism nodes in red (antigen loss, trogocytosis, etc.)
- Edges connecting each antigen to its associated diseases, products, and resistance mechanisms
- A stats line: "Knowledge Graph: 25 targets, 8 toxicities, 10 manufacturing processes, 15+ biomarkers, 6 regulatory products"

**Select:** Entity Type = `Toxicities` to show a different graph view with toxicity nodes linked to biomarkers and management drugs.

**Show:** Scroll down to the **"Cross-Collection Entity Search"** section below the graph.

**Type:** `CD19` in the entity search input (placeholder: "e.g., Yescarta, CD19, FMC63").

**Expected result:** Results grouped by collection — Literature, Trials, Constructs, Assays, etc. — each showing the top hits with relevance scores and text snippets.

**Talking points:**

- "The knowledge graph contains 25 CAR-T target antigens — from well-established CD19 and BCMA to emerging targets like GPRC5D and Claudin18.2."
- "All 6 FDA-approved products are modeled with their complete specifications: target antigen, scFv origin, costimulatory domain, vector type, and approval timeline."
- "8 toxicity profiles cover CRS, ICANS, B-cell aplasia, HLH/MAS, cytopenias, TLS, GvHD, and on-target/off-tumor toxicity."
- "The entity search finds all evidence related to any concept across all 10 collections — a cross-functional lookup."
- "This structured knowledge is automatically injected into every query prompt — Claude knows the entire landscape."

---

### Step 6: Deep Research Mode (2 minutes)

**Click:** The **"Chat"** tab to return to the chat interface.

**Click:** The **"Deep Research Mode"** toggle in the sidebar to enable it. Confirm the badge changes from "QUICK RAG" (blue) to "DEEP RESEARCH" (purple).

**Type this query:**

> What are the key challenges in developing solid tumor CAR-T therapies and what strategies are being explored?

**Expected result:** The UI shows the agent reasoning pipeline in real-time via the status indicator:

1. **"Deep Research: planning search strategy..."** — the agent identifies this as a broad, multi-faceted question. The status area shows:
   - **Strategy:** the search approach the agent chose
   - **Targets:** relevant target antigens identified
   - **Stages:** development stages to search
   - **Sub-questions:** count of decomposed sub-queries
2. **"Deep Research: retrieving evidence..."** — each sub-query gets its own retrieval pass across 11 collections.
3. **Evidence quality** evaluation — the agent checks if evidence is sufficient. If insufficient, it automatically runs supplementary searches with sub-questions and reports "Augmented to: N hits".
4. **"Generating response..."** — Claude produces a comprehensive, multi-section answer streamed token-by-token with citations from every sub-query.

**Talking points:**

- "Deep Research Mode activates the full agent reasoning pipeline — Plan, Decompose, Search, Evaluate, Synthesize."
- "Watch the step-by-step progress in real-time — the agent decides what to search for and how to structure the answer."
- "If the initial search returns insufficient evidence, the agent autonomously expands its search using decomposed sub-questions."
- "This is autonomous reasoning, not a static prompt template."

**Click:** The **"Deep Research Mode"** toggle to disable it (return to Quick RAG mode for the remaining steps).

---

### Step 7: Report Export (1 minute)

**Show:** The export buttons that appear below the most recent chat response.

**Click:** The **"Download PDF"** button.

**Expected result:** The browser downloads a professionally styled PDF report with NVIDIA theming via ReportLab.

**Click:** The **"Download Markdown"** button.

**Expected result:** A `.md` file downloads with the full query, response, evidence sources, and citations.

**Click:** The **"Download JSON"** button.

**Expected result:** A `.json` file downloads with structured data including query metadata, response text, evidence hits, and filter settings.

**Talking points:**

- "Three export formats: Markdown for sharing, JSON for programmatic consumption, PDF for clinical documentation."
- "The PDF is NVIDIA-themed with professional styling via ReportLab."
- "Every report includes: query metadata, RAG-synthesized analysis, evidence sources by collection, knowledge graph context, and citation links."
- "Download the PDF — open it and you'll see a publication-ready report with clickable PubMed links."
- "All three buttons appear after every query — export any answer at any time."

---

### Step 8: Closing Route A (2 minutes)

**Talking points:**

- "10 owned collections, 6,266+ curated vectors, 1 shared genomic evidence collection with 3.56 million vectors."
- "25 target antigens, 6 FDA-approved products, 8 toxicity profiles — all in the knowledge graph."
- "12-16ms retrieval across 11 collections. Full RAG with Claude synthesis in ~24 seconds."
- "Automatic comparative analysis — just say 'compare' and the engine does the rest."
- "Deep Research Mode for complex questions — autonomous agent decomposition and evidence evaluation."
- "This is a complete CAR-T cell therapy intelligence platform — from target discovery through manufacturing to clinical safety."
- "Every answer is grounded in published evidence with clickable citations. No hallucination."
- "And the entire demo was done through a Streamlit interface — no command line, no API calls, no JSON payloads."

---

## Route B: Cross-Platform Integration Demo (25 minutes)

> **Prerequisite:** Complete Route A first, or start from a fresh session with all services running.
>
> This route demonstrates how the CAR-T Intelligence Agent connects to the full HCLS AI Factory — bridging cell therapy intelligence with genomic variant analysis and AI-driven drug discovery.

| Step | UI | Action | Time |
|------|-----|--------|------|
| 1 — Platform Overview | Landing page `:8080` | Show full platform health, architecture | 2 min |
| 2 — Shared Data Layer | Chat tab sidebar `:8521` | Show 10+1 collections, explain shared genomic_evidence | 3 min |
| 3 — Patient Variant Context | Chat tab `:8521` | Query bridging genomic variants and CAR-T therapy | 4 min |
| 4 — CAR-T + Genomics Bridge | Chat tab `:8521` | Query about VCP target and CAR-T vs small molecules | 4 min |
| 5 — Cross-Functional Compare | Chat tab `:8521` | Compare CD19 vs BCMA manufacturing and outcomes | 3 min |
| 6 — Bridge to Stage 2 | RAG Chat `:8501` | Query genomic evidence for antigen loss prediction | 3 min |
| 7 — Bridge to Stage 3 | Drug Discovery UI `:8505` | Show target-to-molecule pipeline | 3 min |
| 8 — The Complete Loop | -- | Walk through full precision medicine pipeline | 3 min |

---

### Step 1: Platform Overview (2 minutes)

**Show:** HCLS AI Factory landing page at http://localhost:8080

**Talking points:**

- "The HCLS AI Factory is a 3-stage precision medicine platform: GPU-accelerated genomics, RAG-grounded target identification, and AI-driven drug discovery."
- "The CAR-T Intelligence Agent extends this platform with specialized cell therapy knowledge."
- "All agents share the same Milvus vector database, BGE embeddings, and Claude LLM — creating a unified intelligence layer."

**Architecture overview:**

```
Stage 1: GPU Genomics (Parabricks)
    | 11.7M variants -> 3.56M quality-filtered
Stage 2: RAG Target Identification (Claude + Milvus)
    | genomic_evidence: 3.56M vectors (shared)
CAR-T Intelligence Agent (11 collections)
    | Target validation + construct intelligence
Stage 3: Drug Discovery (BioNeMo)
    | 100 ranked drug candidates
```

---

### Step 2: Shared Genomic Data Layer (3 minutes)

**Show:** Chat tab sidebar at http://localhost:8521

**Click:** Scroll through the sidebar **Collections** section. Point out all 10 collection checkboxes with their live record counts:

| Collection | Checkbox Label | Vectors |
|---|---|---|
| cart_literature | Literature (5,047) | 5,047 |
| cart_trials | Clinical Trials (973) | 973 |
| cart_constructs | CAR Constructs (6) | 6 |
| cart_assays | Assay Data (45) | 45 |
| cart_manufacturing | Manufacturing (30) | 30 |
| cart_safety | Safety (40) | 40 |
| cart_biomarkers | Biomarkers (43) | 43 |
| cart_regulatory | Regulatory (25) | 25 |
| cart_sequences | Sequences (27) | 27 |
| cart_realworld | Real-World Evidence (30) | 30 |
| genomic_evidence | Genomic Evidence (3,561,170) | 3,561,170 |

**Show:** The **Total** line at the bottom of the collections: "Total: 3,567,436 vectors across 11 collections"

**Talking points:**

- "The CAR-T agent has 10 specialized collections — 6,266 curated vectors covering literature, trials, constructs, assays, manufacturing, safety, biomarkers, regulatory, sequences, and real-world evidence."
- "It also reads from `genomic_evidence` — 3,561,170 vectors from the genomics pipeline. These are the same variant annotations Stage 2 uses."
- "This is the key integration point. The same patient variants that Stage 2 analyzes for drug targets are available to inform CAR-T therapy decisions."

---

### Step 3: Patient Variant Context (4 minutes)

**Show:** Chat tab at http://localhost:8521

**Select:** Verify all 11 collection checkboxes are checked (including Genomic Evidence).

**Type this query:**

> What genomic variants in this patient affect CAR-T therapy response and toxicity prediction?

**Expected result:** The status indicator shows retrieval across multiple collections. Claude searches both the CAR-T therapy collections and the genomic evidence collection. The response includes:

- Evidence cards from `cart_safety` (red badge), `cart_biomarkers` (teal badge), and `genomic_evidence` (cyan badge) with relevance scores
- Cross-functional synthesis connecting variants to CAR-T outcomes

**Talking points:**

- "This query bridges genomic variants and CAR-T intelligence. Claude searches both the therapy collections and the genomic evidence."
- "It identifies variants relevant to CRS prediction — ferritin, IL-6 pathway polymorphisms."
- "It flags variants that could affect T-cell expansion and persistence — immune checkpoint gene variants."
- "This is precision cell therapy — matching the patient's genome to the optimal CAR-T strategy."

---

### Step 4: CAR-T + Genomics Bridge (4 minutes)

**Show:** Chat tab at http://localhost:8521

**Type this query:**

> If the genomics pipeline identified VCP as a drug target for frontotemporal dementia, how would CAR-T approaches compare to small molecule inhibitors for this target?

**Expected result:** Claude pulls from literature on CAR-T for neurodegeneration, constructs targeting intracellular proteins, and manufacturing considerations. Evidence cards show hits from `cart_literature` (blue badge), `cart_constructs` (green badge), and `genomic_evidence` (cyan badge).

**Talking points:**

- "VCP is the demo target from Stage 2 — identified by the genomics pipeline as a frontotemporal dementia drug target."
- "The CAR-T agent evaluates whether cell therapy approaches could work alongside small molecule inhibitors."
- "It pulls from literature on CAR-T for neurodegeneration, constructs targeting intracellular proteins, and manufacturing considerations."
- "This shows the power of cross-functional intelligence — one agent informing another's decision space."

---

### Step 5: Cross-Functional Comparative Query (3 minutes)

**Show:** Chat tab at http://localhost:8521

**Type this query:**

> Compare manufacturing parameters and clinical outcomes for CD19 vs BCMA CAR-T products

**Expected result:** Claude triggers comparative mode (detecting "Compare") and performs dual retrieval for CD19 and BCMA across manufacturing and clinical collections. The expandable "Comparative Evidence" panel shows:

- **CD19** entity header (blue) with evidence cards
- **"-- VS --"** divider in NVIDIA green
- **BCMA** entity header (purple) with evidence cards
- Claude's synthesized comparison table

**Talking points:**

- "Comparative analysis across manufacturing AND clinical collections simultaneously."
- "CD19 products: 4 approved (Kymriah, Yescarta, Tecartus, Breyanzi) — compare transduction efficiency, expansion time, release criteria."
- "BCMA products: 2 approved (Abecma, Carvykti) — different manufacturing challenges, different patient populations."
- "This cross-functional view doesn't exist in any single publication — it's synthesized from 11 collections."

---

### Step 6: Bridge to Stage 2 — RAG Chat (3 minutes)

**Show:** Stage 2 RAG Chat UI at http://localhost:8501

**Type this query:**

> The CAR-T intelligence analysis identified CD19 antigen loss as a primary resistance mechanism. What genomic variants in this patient could predict antigen loss risk?

**Expected result:** The Stage 2 RAG system searches the same `genomic_evidence` collection (3.56M vectors) and returns variant-level evidence for CD19 locus variants, immune checkpoint polymorphisms, and lineage plasticity markers.

**Talking points:**

- "We've transitioned from the CAR-T agent to the Stage 2 RAG pipeline. Same genomic evidence, different analytical lens."
- "Claude is now reasoning over 3.56 million variant annotations — ClinVar pathogenicity, AlphaMissense predictions, and clinical evidence."
- "The CAR-T agent identified the resistance mechanism. Stage 2 RAG identifies the patient-specific genomic risk factors."

---

### Step 7: Bridge to Stage 3 — Drug Discovery (3 minutes)

**Show:** Drug Discovery UI at http://localhost:8505

**Talking points:**

- "When the CAR-T agent identifies resistance mechanisms — like BCMA downregulation in myeloma — those targets feed into Stage 3."
- "BioNeMo generates small molecule candidates that could be combined with CAR-T therapy."
- "Combination strategies: CAR-T for bulk disease + small molecule for resistant clones."
- "The platform enables multi-modal therapeutic design — not just cell therapy OR drug discovery, but both together."

---

### Step 8: The Complete Loop (3 minutes)

**Talking points:**

Walk through the full precision medicine loop:

```
Patient DNA (Illumina Sequencing)
    |
Stage 1: GPU Genomics (Parabricks)
    | 11.7M variants called
Stage 2: RAG Target Identification (Claude + Milvus)
    | VCP identified as drug target
    | genomic_evidence: 3.56M vectors (shared)
    |---> CAR-T Intelligence Agent
    |     - Evaluates cell therapy approaches for target
    |     - Cross-references 6,266 therapy-specific vectors
    |     - Identifies optimal construct design + manufacturing
    |---> Imaging Intelligence Agent
    |     - Cross-modal triggers from imaging findings
    |     - Connects phenotype to genotype
    +---> Stage 3: Drug Discovery (BioNeMo)
          - 100 novel drug candidates generated
          - Docking + drug-likeness scoring

Combined Clinical Output:
    -> FHIR R4 Reports (imaging)
    -> PDF Reports (CAR-T + drug discovery)
    -> JSON data for dashboards
```

- "Patient DNA to therapeutic strategy. Genomics, cell therapy intelligence, imaging AI, and drug discovery — all on one platform."
- "Every agent sees the same genomic truth. Every answer is grounded in evidence."
- "This is the future of precision medicine — multi-modal, multi-agent, evidence-grounded."

**Closing Route B:**

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

### Streamlit UI Not Loading

```bash
# Check container status
docker compose ps | grep streamlit

# Check Streamlit logs
docker compose logs cart-streamlit

# Restart the UI container
docker compose restart cart-streamlit
```

Confirm the UI is accessible at http://localhost:8521. If port 8521 is occupied, check for conflicting services.

### PDF Export Issues

PDF generation requires ReportLab. If PDF export fails:

```bash
pip install reportlab
```

---

## Quick Reference

| Resource | URL |
|---|---|
| CAR-T Agent UI | http://localhost:8521 |
| CAR-T Agent API docs | http://localhost:8522/docs |
| RAG Chat (Stage 2) | http://localhost:8501 |
| Drug Discovery UI | http://localhost:8505 |
| Landing Page | http://localhost:8080 |
| Grafana Monitoring | http://localhost:3000 |
| Milvus | http://localhost:19530 |
| Attu (Milvus UI) | http://localhost:8000 |

---

## Appendix: API Reference (Developer Use)

> All CAR-T Intelligence Agent endpoints are documented at http://localhost:8522/docs (FastAPI auto-generated Swagger UI). The following curl commands are provided for scripting, automated testing, and integration development. These are **not** part of the live demo — all demo interaction uses the Streamlit UI.

### Health and Status

```bash
# Health check — returns collection counts and total vectors
curl -s http://localhost:8522/health | python3 -m json.tool

# Knowledge graph statistics
curl -s http://localhost:8522/knowledge/stats | python3 -m json.tool

# List all collections with vector counts
curl -s http://localhost:8522/collections | python3 -m json.tool

# Prometheus metrics
curl -s http://localhost:8522/metrics
```

### RAG Query

```bash
# Full RAG query with Claude synthesis
curl -s -X POST http://localhost:8522/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Why do CD19 CAR-T therapies fail in relapsed B-ALL?",
    "target_gene": "CD19"
  }' | python3 -m json.tool
```

Response fields:

| Field | Description |
|---|---|
| `answer` | Claude-synthesized response with citations |
| `sources` | Array of evidence hits with collection, score, and text snippet |
| `follow_up_questions` | Suggested next queries |
| `confidence` | Overall confidence score (0-1) |
| `processing_time_ms` | End-to-end latency |

### Full Query (with comparative detection)

```bash
# Triggers automatic comparative analysis when "compare", "vs", "versus" detected
curl -s -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Compare 4-1BB vs CD28 costimulatory domains for DLBCL"
  }' | python3 -m json.tool
```

### Evidence-Only Search (no LLM synthesis)

```bash
# Raw vector search — returns ranked evidence without Claude synthesis
curl -s -X POST http://localhost:8522/search \
  -H "Content-Type: application/json" \
  -d '{
    "question": "BCMA CAR-T resistance mechanisms in multiple myeloma",
    "target_antigen": "BCMA",
    "top_k": 5
  }' | python3 -m json.tool
```

### Find Related Entities

```bash
curl -s -X POST http://localhost:8522/find-related \
  -H "Content-Type: application/json" \
  -d '{
    "entity": "CD19",
    "entity_type": "gene"
  }' | python3 -m json.tool
```

### Report Generation

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

# PDF report (binary download)
curl -s -X POST http://localhost:8522/api/report/pdf \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What manufacturing parameters predict CAR-T clinical response?"
  }' --output cart_manufacturing_report.pdf
```

### Performance Benchmarks

| Operation | Typical Latency |
|---|---|
| Embedding (BGE-small-en-v1.5) | < 5 ms |
| 11-collection parallel search | 12-16 ms |
| Query expansion + re-search | 8-12 ms |
| Merge + deduplicate + rank | < 1 ms |
| **Total retrieval (no LLM)** | **~25 ms** |
| Full RAG with Claude synthesis | ~24 sec |

---

*HCLS AI Factory — Apache 2.0 | February 2026*
