# CAR-T Intelligence Agent -- Demo Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 2.0.0
**License:** Apache 2.0

**Estimated Demo Time:** 30 minutes (expandable to 45 with full API walkthrough)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Pre-Demo Checklist](#2-pre-demo-checklist)
3. [Opening Hook (~2 min)](#3-opening-hook-2-min)
4. [Scenario 1: Basic Query (~5 min)](#4-scenario-1-basic-query-5-min)
5. [Scenario 2: Comparative Analysis (~5 min)](#5-scenario-2-comparative-analysis-5-min)
6. [Scenario 3: Cross-Functional Query (~5 min)](#6-scenario-3-cross-functional-query-5-min)
7. [Scenario 4: Safety Intelligence (~3 min)](#7-scenario-4-safety-intelligence-3-min)
8. [Scenario 5: Genomic Bridge (~3 min)](#8-scenario-5-genomic-bridge-3-min)
9. [Scenario 6: Patent and IP (~2 min)](#9-scenario-6-patent-and-ip-2-min)
10. [Scenario 7: Immunogenicity (~2 min)](#10-scenario-7-immunogenicity-2-min)
11. [Advanced Features (~5 min)](#11-advanced-features-5-min)
12. [Closing: The Big Picture](#12-closing-the-big-picture)
13. [Troubleshooting](#13-troubleshooting)
14. [Quick Reference Card](#14-quick-reference-card)

---

## 1. Overview

### What This Demo Shows

The CAR-T Intelligence Agent is a cross-functional AI system that searches **3,567,622 vectors across 11 specialized data collections** to answer questions about CAR-T cell therapy development. It spans the entire lifecycle from patient genomics and target identification through manufacturing, clinical trials, post-market safety, and real-world outcomes.

This is not a generic chatbot. It is a domain-specific intelligence platform that:

- Searches 11 Milvus collections **in parallel** with a single query
- Augments responses with a **structured knowledge graph** (34 target antigens, 17 toxicity profiles, 20 manufacturing processes, 23 biomarkers, 6 regulatory products, 6 immunogenicity topics)
- Produces **clickable PubMed and ClinicalTrials.gov citations**
- Auto-detects comparative queries and generates **structured side-by-side analysis**
- Bridges to **3.5 million genomic variants** from the HCLS AI Factory genomics pipeline
- Exports results to **Markdown, JSON, and PDF**

### Who This Demo Is For

- Pharmaceutical and biotech executives evaluating AI platforms
- Cell therapy researchers and development teams
- Regulatory affairs and pharmacovigilance professionals
- Clinical operations and manufacturing leaders
- IT and data science teams evaluating RAG architectures
- Conference and trade show audiences (GTC, ASGCT, ASH, ASCO)

### Key Message

> "What if every question about your CAR-T program could instantly search across published literature, clinical trials, manufacturing records, safety data, biomarkers, regulatory filings, molecular sequences, real-world evidence, and 3.5 million patient genomic variants -- all at once?"

---

## 2. Pre-Demo Checklist

Complete these steps **at least 15 minutes before the demo**.

### 2.1 Verify Milvus Is Running

```bash
curl -s http://localhost:19530/v1/vector/collections | python3 -m json.tool | head -5
```

**Expected:** JSON response listing collections. If Milvus is not running, start it:

```bash
# If using Docker
docker start milvus-standalone

# Or if using the HCLS AI Factory docker-compose
cd /home/adam/projects/hcls-ai-factory && docker compose up -d milvus
```

### 2.2 Verify the ANTHROPIC_API_KEY Is Available

```bash
# Check environment variable
echo $ANTHROPIC_API_KEY | head -c 10

# Or check the .env file
head -1 /home/adam/projects/hcls-ai-factory/rag-chat-pipeline/.env
```

**Expected:** The key should start with `sk-ant-`. If missing, the agent will run in scaffold mode (retrieval works but LLM synthesis will not).

### 2.3 Start the Streamlit UI

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/cart_intelligence_agent
streamlit run app/cart_ui.py --server.port 8521
```

**Expected:** Terminal shows "You can now view your Streamlit app in your browser" with the URL `http://localhost:8521`.

### 2.4 Start the FastAPI Server (Optional, for API Demo)

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/cart_intelligence_agent
uvicorn api.main:app --host 0.0.0.0 --port 8522 --reload
```

**Expected:** Terminal shows "Uvicorn running on http://0.0.0.0:8522".

### 2.5 Verify System Health

Open the browser to `http://localhost:8521`. The sidebar should show all 11 collections with record counts:

| Collection | Expected Count |
|---|---|
| Literature | 5,047 |
| Clinical Trials | 973 |
| CAR Constructs | 41 |
| Assay Data | 75 |
| Manufacturing | 56 |
| Safety | 71 |
| Biomarkers | 60 |
| Regulatory | 40 |
| Sequences | 40 |
| Real-World Evidence | 54 |
| Genomic Evidence | 3,561,170 |
| **Total** | **3,567,622** |

If any collection shows 0, run the appropriate seed script from the `scripts/` directory.

### 2.6 Pre-Warm the Embedding Model

The first query takes a few seconds longer because the BGE-small-en-v1.5 model needs to load into memory. Submit one throwaway query before the audience arrives:

> Type: "test query" and press Enter. Wait for the response, then clear the chat.

### 2.7 Browser Setup

- Open `http://localhost:8521` in a clean browser tab
- Ensure the sidebar is expanded (it should be by default)
- Set zoom to 100% or 110% for readability on projectors
- If using a second screen, open `http://localhost:8522/docs` for the FastAPI Swagger UI

---

## 3. Opening Hook (~2 min)

### The Problem Statement

> **Say this:**
> "CAR-T cell therapy is one of the most complex therapeutic modalities ever developed. A single CAR-T program generates data across at least a dozen different domains -- target biology, construct design, manufacturing, preclinical testing, clinical trials, safety monitoring, biomarkers, regulatory filings, molecular sequences, real-world outcomes, and patient genomics. Today, this data lives in silos. A manufacturing scientist cannot easily search clinical trial results. A regulatory affairs specialist cannot easily correlate safety signals with molecular design choices. A clinical researcher cannot easily link a patient's genomic profile to predicted CAR-T response."

> **Then say:**
> "What if you could break down all of those silos with a single question?"

### Show the Numbers

Point to the sidebar. Read the total aloud:

> "This system has indexed **3,567,622 vectors** across **11 specialized collections**. That includes over 5,000 published research papers, nearly 1,000 clinical trials, safety reports, biomarker data, manufacturing records, molecular sequences, real-world registry outcomes, and -- here is the bridge to precision medicine -- **3.5 million genomic variants** from actual patient sequencing data processed through our Parabricks genomics pipeline."

### Show the Collection List

Scroll down the sidebar checkboxes and briefly name each:

> "Literature, Clinical Trials, CAR Constructs, Assays, Manufacturing, Safety, Biomarkers, Regulatory, Sequences, Real-World Evidence, and Genomic Evidence. Every single one of these is searched simultaneously with every query."

---

## 4. Scenario 1: Basic Query (~5 min)

### Goal

Demonstrate the core RAG pipeline: multi-collection search, knowledge graph augmentation, streaming LLM response, evidence panel with citations.

### Query to Type

```
Why do CD19 CAR-T therapies fail in relapsed B-ALL?
```

**Alternative (from the demo buttons):** Click the first demo query button in the sidebar: "Why do CD19 CAR-T therapies fail in relapsed B-ALL?"

### What to Point Out During the Search

While the status indicator shows "Searching across CAR-T data sources...":

> "Watch the search status. It is searching all 11 collections in parallel using a ThreadPoolExecutor. The query is embedded with BGE-small-en-v1.5 -- a 384-dimensional vector -- and compared against 3.5 million vectors using cosine similarity. This takes about 20-30 milliseconds for the retrieval step."

### What to Point Out in the Results

1. **Search metrics:** When the status expands, read the counts:
   > "It found [X] results across [Y] collections in [Z] milliseconds."

2. **Collection distribution:** Note which collections contributed:
   > "Look at the breakdown -- we have hits from Literature, Trials, Constructs, Assays, Safety, and Biomarkers. This is a cross-functional answer, not just a literature search."

3. **Streaming response:** As the answer streams in, point out:
   > "The response is streaming from Claude Sonnet 4.6. Notice the citations -- those are clickable links. Every claim is grounded in the evidence retrieved from our vector database."

4. **Click a PubMed link:** Click one of the `[Literature:PMID XXXXXXXX]` links to show it opens the actual PubMed abstract in a new tab.

5. **Click a ClinicalTrials.gov link:** Click one of the `[Trial:NCTXXXXXXXX]` links to show it opens the actual trial record.

6. **Evidence panel:** Expand the "Evidence Sources" expander below the response:
   > "Here is every piece of evidence the system retrieved. Each card shows the collection it came from -- color-coded by type -- the cosine similarity score, the relevance badge, and a snippet of the actual text."

### Talking Points

- **Knowledge graph augmentation:** "The system detected 'CD19' in the query and automatically injected structured knowledge about CD19: known resistance mechanisms (antigen loss, lineage switch, trogocytosis, alternative splicing), approved products (Kymriah, Yescarta, Tecartus, Breyanzi), and the toxicity profile (CRS 30-90%, ICANS 20-65%). This structured data complements the retrieved evidence."

- **Query expansion:** "Behind the scenes, the query expansion engine detected 'CD19' and expanded the search to include related terms like B-ALL, DLBCL, tisagenlecleucel, axicabtagene ciloleucel, FMC63. This improves recall across all collections."

- **Relevance scoring:** "Notice the green, yellow, and grey badges on the evidence cards. Green means high relevance (cosine similarity 0.75 or above), yellow is medium (0.60-0.75), and grey is lower relevance. The system prioritizes high-relevance citations in the generated response."

---

## 5. Scenario 2: Comparative Analysis (~5 min)

### Goal

Demonstrate auto-detected comparative analysis with structured side-by-side output and dual retrieval.

### Query to Type

```
Compare 4-1BB vs CD28 costimulatory domains
```

### What to Point Out During the Search

> "The system auto-detected this as a comparative query because of the word 'compare' and the 'vs' keyword. It parsed two entities: 4-1BB and CD28. Each entity is resolved against our knowledge graph -- 4-1BB maps to CD137/TNFRSF9, and CD28 is a costimulatory domain used in Yescarta and Tecartus. The system then runs **two separate retrieval pipelines**, one focused on each entity."

Watch for the status message that says: "Comparative analysis: **4-1BB (CD137)** vs **CD28**"

### What to Point Out in the Results

1. **Structured comparison output:** The response should include:
   - A **comparison table** in markdown format with dimensions as rows and the two domains as columns
   - **Advantages** of each domain (bulleted)
   - **Limitations** of each domain (bulleted)
   - A **clinical context** paragraph

   > "This is not a free-form essay. The system instructed Claude to produce a structured comparison with a table, advantages, limitations, and clinical context. Every claim is cited."

2. **Comparative evidence panel:** Expand the evidence panel. It should show:
   - A **blue header** for 4-1BB evidence
   - A **"-- VS --"** divider in green
   - A **purple header** for CD28 evidence

   > "The evidence is grouped by entity. You can see which papers and trials support claims about each costimulatory domain independently."

3. **Products mentioned:** The response should reference specific FDA products:
   > "Notice it mentions Kymriah and Breyanzi (4-1BB products) versus Yescarta and Tecartus (CD28 products), with specific trial names like ELIANA, ZUMA-1, and TRANSCEND."

### Talking Points

- **Entity resolution:** "The system resolves 'Kymriah' to CD19, '4-1BB' to a costimulatory domain, 'BCMA' to a target antigen. It understands the CAR-T taxonomy because of the knowledge graph with 34 antigens, 6 products, and 54 aliases."

- **Other comparisons you can try:**
  - "Kymriah versus Yescarta" (product-level comparison, resolves to CD19 with different constructs)
  - "Compare CRS and ICANS" (toxicity profile comparison)
  - "Compare lentiviral vs retroviral transduction" (manufacturing comparison)

---

## 6. Scenario 3: Cross-Functional Query (~5 min)

### Goal

Demonstrate the system's ability to span multiple development stages in a single query, connecting manufacturing to clinical outcomes. Show knowledge graph augmentation in depth.

### Query to Type

```
What manufacturing parameters predict clinical response?
```

### What to Point Out During the Search

> "This is a deliberately broad query. It does not name a specific target antigen or product. The system needs to find evidence spanning manufacturing processes, clinical trial outcomes, and biomarker data to answer this question. Watch which collections contribute."

### What to Point Out in the Results

1. **Multi-domain evidence:** The search should return hits from at least 4-5 collections:
   - **Manufacturing:** T-cell expansion parameters, transduction efficiency, VCN
   - **Literature:** Correlative studies linking product attributes to outcomes
   - **Biomarkers:** Tcm percentage, CD4:CD8 ratio, exhaustion markers
   - **Trials:** Responder versus non-responder analyses from ELIANA, ZUMA-1, KarMMa
   - **Constructs:** Product-specific manufacturing differences

   > "This is why we built 11 collections instead of one. A single literature collection would give you published abstracts. But by having separate manufacturing and biomarker collections, the system retrieves specific process parameters alongside the clinical evidence that validates them."

2. **Knowledge graph contribution:** The system should inject manufacturing knowledge:
   > "The knowledge graph detected 'manufacturing' in the query and injected structured data about lentiviral transduction parameters, expansion protocols, and release testing criteria. This gives Claude specific parameter names and thresholds to reference."

3. **Specific parameters to listen for:**
   - T-cell fitness / central memory percentage (Tcm >40%)
   - CD4:CD8 ratio (defined composition in Breyanzi)
   - Vector copy number (VCN <5 copies/cell)
   - Expansion fold and duration (rapid 6-day Kite process vs. 9-12 day standard)
   - Post-thaw viability (>70%)
   - Vein-to-vein time (3-6 weeks centralized vs. 3-7 days point-of-care)

### Talking Points

- **Breaking silos:** "In most organizations, the manufacturing team has their data in batch records, the clinical team has their data in EDC systems, and the correlative science team has their data in lab notebooks. This system unifies all of those data types into a single searchable intelligence layer."

- **Actionable insights:** "The answer is not just 'manufacturing matters.' It tells you which specific parameters to optimize: T-cell fitness, Tcm frequency, and the ratio of CD4 to CD8 cells in the final product. These are actionable insights that a manufacturing scientist can use today."

---

## 7. Scenario 4: Safety Intelligence (~3 min)

### Goal

Demonstrate safety and biomarker collection retrieval, the target antigen filter, and clinical grading system knowledge.

### Query to Type

```
Which biomarkers best predict CRS severity?
```

### What to Point Out

1. **Safety + Biomarker collection hits:** The evidence panel should show cards from both the Safety collection (red badges) and the Biomarker collection (teal badges):
   > "We have dedicated Safety and Biomarker collections. The Safety collection contains pharmacovigilance data with grading information. The Biomarker collection contains predictive and pharmacodynamic markers with clinical cutoffs and assay methods."

2. **Specific biomarkers:** The response should discuss:
   - Serum ferritin (>500 mg/L pre-infusion predicts grade 3+ CRS)
   - C-reactive protein (>200 mg/L within 72 hours)
   - IL-6 (>1000 pg/mL, target of tocilizumab)
   - Soluble IL-2 receptor (sCD25)
   - Peak CAR-T expansion (Cmax)

3. **Toxicity management:** The knowledge graph injects CRS management protocols:
   > "The system injected the full CRS grading system (Lee 2014 / ASTCT 2019) and management ladder: tocilizumab first-line for grade 2+, corticosteroids for grade 3+, siltuximab for refractory cases, and emerging agents like anakinra."

### Demonstrate the Target Antigen Filter

While the response is visible, go to the sidebar and change the **Target Antigen Filter** dropdown from "All Targets" to "CD19". Then ask a follow-up:

```
What are the long-term safety signals for CD19 CAR-T products?
```

> "Now the search is filtered to CD19-specific evidence only. This narrows the results to Kymriah, Yescarta, Tecartus, and Breyanzi safety data, excluding BCMA products."

---

## 8. Scenario 5: Genomic Bridge (~3 min)

### Goal

Demonstrate the unique integration with the 3.5 million genomic variants from the HCLS AI Factory genomics pipeline.

### Query to Type

```
What genomic variants in CD19 or BCMA pathway genes affect CAR-T response?
```

### What to Point Out

1. **Genomic Evidence collection:** The evidence panel should show hits from the Genomic Evidence collection (purple badges) alongside Literature and other collections:
   > "Here is where the CAR-T Intelligence Agent connects to the broader HCLS AI Factory platform. Those 3.5 million vectors in the genomic_evidence collection come from actual patient VCF data processed through our Parabricks genomics pipeline. They include ClinVar clinical significance annotations and AlphaMissense pathogenicity predictions."

2. **Variant-level detail:** Genomic evidence cards show:
   - Gene name
   - Consequence type (missense, frameshift, etc.)
   - Impact level (HIGH, MODERATE, LOW, MODIFIER)
   - Clinical significance from ClinVar
   - AlphaMissense pathogenicity class

3. **Cross-domain synthesis:** The LLM response should connect genomic variants to CAR-T therapy mechanisms:
   > "The system is connecting genetic variants in CD19 pathway genes to known resistance mechanisms. For example, CD19 mutations can cause antigen loss, which is the primary mechanism of relapse in 20-30% of patients treated with CD19-directed CAR-T therapy."

### Talking Points

- **Precision medicine vision:** "This is the precision medicine vision: patient DNA to drug candidates. The genomics pipeline identifies variants. The RAG pipeline provides clinical context. And now the CAR-T agent connects those variants to specific cell therapy design and treatment decisions."

- **Pipeline integration:** "The genomic_evidence collection was created by Stage 2 of the HCLS AI Factory pipeline. The CAR-T agent reads it as a shared resource. No data was duplicated -- it is the same Milvus instance, the same embedding model, the same vectors."

---

## 9. Scenario 6: Patent and IP (~2 min)

### Goal

Demonstrate patent-related search capability from the literature collection.

### Query to Type

```
What patents cover bispecific CAR-T constructs targeting CD19 and CD22?
```

### What to Point Out

1. **Patent literature:** The Literature collection includes patent records alongside PubMed abstracts. Evidence cards with patent IDs demonstrate the breadth of the knowledge base.

2. **Bispecific construct knowledge:** The query expansion engine detects "bispecific" and expands to: tandem CAR, dual-targeting, bivalent CAR, OR-gate logic, CD19/CD22, bicistronic, loop CAR, split CAR.

3. **IP landscape awareness:** The response should discuss:
   - Dual-targeting strategies to prevent antigen escape
   - Tandem CAR designs (two scFvs in a single chain)
   - Bicistronic approaches (two separate CARs from one vector)
   - Key academic and commercial entities with IP in this space

### Talking Points

- **Competitive intelligence:** "This is useful not just for scientific questions but for IP and competitive landscape analysis. A single query can surface patent filings, published research, clinical trial registrations, and construct designs related to a specific technology."

---

## 10. Scenario 7: Immunogenicity (~2 min)

### Goal

Demonstrate the immunogenicity and HLA knowledge domain, including humanization strategies and anti-drug antibody (ADA) risk.

### Query to Type

```
How does scFv humanization reduce immunogenicity risk in CAR-T therapy?
```

### What to Point Out

1. **Immunogenicity knowledge graph:** The system injects detailed knowledge about:
   - Murine scFv immunogenicity (FMC63 in Kymriah/Yescarta, ADA incidence 3-8%)
   - Humanization strategies (CDR grafting, framework shuffling, deimmunization)
   - HLA-restricted T-cell epitopes (HLA-DRB1*04:01 and HLA-DRB1*15:01 as high-risk alleles)
   - Computational tools (NetMHCIIpan, EpiMatrix, IEDB)

2. **Sequence collection hits:** Evidence cards from the Sequence collection (indigo badges) show molecular data: scFv clone names, binding affinities, species of origin, and immunogenicity risk ratings.

3. **Biomarker connection:** Biomarker collection hits for ADA monitoring complement the molecular design data.

### Talking Points

- **Translational depth:** "This query spans three collections: Sequences for molecular design data, Biomarkers for ADA monitoring, and Literature for published humanization studies. The knowledge graph adds specific HLA alleles, prediction tools, and ADA incidence rates. This is the kind of cross-functional insight that would normally require consulting three different teams."

---

## 11. Advanced Features (~5 min)

### 11.1 Collection Filters

**Demonstrate:** In the sidebar, uncheck all collections except "Clinical Trials" and "Safety." Then ask:

```
What are the key safety findings from pivotal CAR-T trials?
```

> "By toggling collection filters, you can focus the search on specific data types. This is useful when a regulatory reviewer wants only trial and safety data without the noise of preclinical literature."

Re-enable all collections when done.

### 11.2 Date Range Filter

**Demonstrate:** Check "Apply date filter" in the sidebar. Set the range to 2022-2026. Then ask:

```
What are the latest advances in CAR-T manufacturing?
```

> "The date filter restricts results to publications and trials from the last few years, surfacing only the most recent advances."

Uncheck "Apply date filter" when done.

### 11.3 Development Stage Filter

**Demonstrate:** Change the "Development Stage" dropdown to "Clinical." Then ask:

```
How does T-cell exhaustion affect persistence?
```

> "Filtering by development stage focuses the search on clinical-stage evidence, filtering out preclinical or manufacturing-focused records."

Reset to "All Stages" when done.

### 11.4 Export Results

After any query with results, point out the three download buttons below the response:

1. Click **"Download Markdown"** -- opens/saves a `.md` file with the full report including evidence tables
2. Click **"Download JSON"** -- opens/saves a `.json` file with structured data suitable for programmatic consumption
3. Click **"Download PDF"** -- opens/saves a styled PDF with NVIDIA-themed formatting, evidence tables, and clickable citation links

> "Every query result can be exported in three formats. The PDF is presentation-ready with NVIDIA green branding. The JSON is machine-readable for integration with other systems. The Markdown is ideal for documentation or sharing in collaborative tools."

### 11.5 Conversation Memory (Multi-Turn)

**Demonstrate:** After the response to a previous query is visible, type a follow-up question:

```
What about in the pediatric population specifically?
```

> "The system maintains conversation memory across turns. It injected the context from the previous exchange into the current query so Claude understands 'what about' refers to the previous topic. This enables natural multi-turn conversations without repeating context."

### 11.6 Deep Research Mode

**Demonstrate:** Toggle "Deep Research Mode" ON in the sidebar. Then ask:

```
Why do CD19 CAR-T therapies fail in relapsed B-ALL?
```

> "Deep Research mode activates the autonomous agent pipeline. Instead of a single retrieval, the agent plans a search strategy, evaluates evidence quality, and decomposes complex questions into sub-queries for additional coverage. Watch the status messages -- it shows the strategy, identified targets, development stages, sub-questions, and evidence quality assessment."

Toggle Deep Research Mode OFF when done.

### 11.7 Knowledge Graph Tab

Click the **"Knowledge Graph"** tab at the top of the page.

> "This is an interactive visualization of the CAR-T knowledge graph. You can explore target antigens and their relationships to diseases, products, and resistance mechanisms."

**Demonstrate:**
1. Select "Target Antigens" from the entity type dropdown -- shows nodes for CD19, BCMA, etc. with connections to diseases and products
2. Select "Toxicities" -- shows CRS, ICANS with connections to biomarkers and management drugs
3. Select "Regulatory" -- shows FDA-approved products with connections to indications and designations

Also demonstrate the **Cross-Collection Entity Search** at the bottom:

Type "Yescarta" in the entity search box:

> "This searches for everything related to Yescarta across all 11 collections -- literature, trials, constructs, safety records, regulatory milestones, manufacturing data. It is an entity-centric view of the knowledge base."

### 11.8 Image Analysis Tab

Click the **"Image Analysis"** tab.

> "You can upload a slide image or document screenshot. The agent uses Claude Vision to extract claims from the image, then searches the knowledge base to verify each claim against the 3.5 million indexed vectors. This is useful for validating slide decks, poster presentations, or regulatory submissions."

### 11.9 API Demo with curl

If the FastAPI server is running on port 8522, demonstrate the REST API:

**Health Check:**

```bash
curl -s http://localhost:8522/health | python3 -m json.tool
```

Expected output:
```json
{
    "status": "healthy",
    "collections": 11,
    "total_vectors": 3567622
}
```

**List Collections:**

```bash
curl -s http://localhost:8522/collections | python3 -m json.tool
```

Expected output: All 11 collections with record counts.

**Evidence-Only Search (Fast, No LLM):**

```bash
curl -s -X POST http://localhost:8522/search \
  -H "Content-Type: application/json" \
  -d '{"question": "CD19 antigen loss resistance mechanism"}' \
  | python3 -m json.tool | head -30
```

> "This is the `/search` endpoint -- evidence retrieval only, no LLM generation. It returns in under 100 milliseconds. Useful for building downstream applications that need fast evidence lookup."

**Full RAG Query (Retrieve + LLM Synthesis):**

```bash
curl -s -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What are the key resistance mechanisms for BCMA CAR-T therapy?"}' \
  | python3 -m json.tool | head -50
```

> "This is the full RAG pipeline via the API: retrieve evidence, augment with the knowledge graph, and synthesize a response with Claude. This takes about 20-25 seconds, dominated by LLM generation."

**Filtered Query (Target Antigen + Year Range):**

```bash
curl -s -X POST http://localhost:8522/search \
  -H "Content-Type: application/json" \
  -d '{
    "question": "CRS management tocilizumab",
    "target_antigen": "CD19",
    "year_min": 2022,
    "year_max": 2026
  }' | python3 -m json.tool | head -30
```

> "The API supports all the same filters as the UI: target antigen, year range, and specific collection selection."

**Cross-Collection Entity Search:**

```bash
curl -s -X POST http://localhost:8522/find-related \
  -H "Content-Type: application/json" \
  -d '{"entity": "Kymriah", "top_k": 3}' \
  | python3 -m json.tool | head -40
```

> "The `find-related` endpoint searches for a single entity across all 11 collections. This powers the entity search in the Knowledge Graph tab."

**Knowledge Graph Statistics:**

```bash
curl -s http://localhost:8522/knowledge/stats | python3 -m json.tool
```

Expected output:
```json
{
    "target_antigens": 34,
    "targets_with_approved_products": 2,
    "toxicity_profiles": 17,
    "manufacturing_processes": 20,
    "biomarkers": 23,
    "regulatory_products": 6,
    "immunogenicity_topics": 6
}
```

**Prometheus Metrics:**

```bash
curl -s http://localhost:8522/metrics
```

> "The `/metrics` endpoint exposes Prometheus-compatible counters for request volume, error rates, and per-collection vector counts. This plugs directly into the Grafana dashboard in the HCLS AI Factory monitoring stack."

---

## 12. Closing: The Big Picture

### How This Fits Into the HCLS AI Factory

> "The CAR-T Intelligence Agent is the fourth stage of the HCLS AI Factory platform. Let me show you the full pipeline."

Draw or show this architecture on a whiteboard or slide:

```
Stage 1: Genomics Pipeline (Parabricks)
  Patient FASTQ --> BWA-MEM2 --> DeepVariant --> 11.7M Variants
  Time: 120-240 min on DGX Spark

Stage 2: RAG/Chat Pipeline (Milvus + Claude)
  VCF --> ClinVar (2.7M) + AlphaMissense (71M) --> 3.5M Vectors in Milvus
  Interactive chat with variant-level precision medicine intelligence

Stage 3: Drug Discovery Pipeline (BioNeMo)
  Target --> MolMIM (molecule generation) --> DiffDock (docking) --> RDKit (scoring)
  Time: 8-16 min for 10 candidate molecules

Stage 4: CAR-T Intelligence Agent   <-- YOU ARE HERE
  11 collections, 3.5M+ vectors, cross-functional CAR-T lifecycle intelligence
  Bridges to genomic evidence from Stage 2
```

### The Hardware Story

> "All of this runs on a single NVIDIA DGX Spark. That is a $4,699 workstation with a GB10 GPU, 128 GB of unified memory, and 20 ARM cores. There is no cloud dependency for the compute -- this is on-premise, sovereign AI for healthcare and life sciences."

### The Open-Source Story

> "The entire platform is Apache 2.0 licensed and open-source. Every line of code, every configuration, every knowledge graph entry is available on GitHub. There are no proprietary locks."

### The Generalizability Story

> "The key architectural insight is that this platform is **not disease-specific**. The same Milvus instance, the same embedding model, the same RAG architecture that powers the VCP/Frontotemporal Dementia drug discovery pipeline now powers a completely different therapeutic modality: CAR-T cell therapy. By changing the knowledge graph, query expansion maps, and collection schemas, you can adapt this to any therapeutic area -- immuno-oncology, gene therapy, rare diseases, neuroscience. The infrastructure is the same."

---

## 13. Troubleshooting

### Problem: "Failed to initialize" Error on UI Load

**Cause:** Milvus is not running or not accessible on port 19530.

**Fix:**
```bash
# Check if Milvus is running
docker ps | grep milvus

# Start Milvus if stopped
docker start milvus-standalone

# Or restart via docker-compose
cd /home/adam/projects/hcls-ai-factory && docker compose up -d milvus

# Wait 10-15 seconds for Milvus to be ready, then refresh the UI
```

### Problem: "LLM generation error" or Scaffold Mode

**Cause:** ANTHROPIC_API_KEY is not set or is invalid.

**Fix:**
```bash
# Set the key directly
export ANTHROPIC_API_KEY="sk-ant-..."

# Or ensure it exists in the .env file
cat /home/adam/projects/hcls-ai-factory/rag-chat-pipeline/.env | grep ANTHROPIC_API_KEY
```

The UI will auto-load the key from `rag-chat-pipeline/.env` if the environment variable is not set.

### Problem: Slow First Query (~10-15 Seconds Before Search Starts)

**Cause:** The BGE-small-en-v1.5 embedding model is loading into memory on the first query. Subsequent queries will be fast.

**Fix:** Pre-warm the model before the demo by submitting a test query (see Pre-Demo Checklist section 2.6).

### Problem: Collection Shows 0 Records in the Sidebar

**Cause:** The collection was not seeded after creation.

**Fix:** Run the appropriate seed script:
```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/cart_intelligence_agent

# Literature
python3 scripts/ingest_pubmed.py --max-results 5000

# Clinical Trials
python3 scripts/ingest_clinical_trials.py --max-results 1500

# Constructs + Assays + Manufacturing + Safety + Biomarkers + Regulatory + Sequences + Real-World
python3 scripts/setup_collections.py --seed-constructs
python3 scripts/seed_assays.py
python3 scripts/seed_manufacturing.py
python3 scripts/seed_safety.py
python3 scripts/seed_biomarkers.py
python3 scripts/seed_regulatory.py
python3 scripts/seed_sequences.py
python3 scripts/seed_realworld.py
```

### Problem: Comparative Query Falls Back to Normal Mode

**Cause:** The entity parser could not resolve one or both entities in the "X vs Y" query.

**Fix:** Use recognized entity names. Supported formats:
- Target antigens: CD19, BCMA, CD22, CD20, etc.
- Products: Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti
- Generic names: tisagenlecleucel, axicabtagene ciloleucel, etc.
- Costimulatory domains: 4-1BB, CD28
- Toxicities: CRS, ICANS

The query must contain "compare," "vs," "versus," or "comparing" for comparative detection.

### Problem: API Server Returns 503

**Cause:** The engine did not initialize correctly during FastAPI startup.

**Fix:**
```bash
# Check the server logs for the specific error
# Usually Milvus connection or embedding model loading issue

# Restart the server
uvicorn api.main:app --host 0.0.0.0 --port 8522 --reload
```

### Problem: Evidence Panel Shows No Results

**Cause:** The query may be too specific or the score threshold is filtering out low-relevance results. The default threshold is 0.4 cosine similarity.

**Fix:** Try a broader query, or check if the relevant collections are enabled in the sidebar.

### Problem: PDF Export Fails

**Cause:** The `reportlab` library is not installed.

**Fix:**
```bash
pip install reportlab
```

### Problem: Knowledge Graph Tab Shows "Install pyvis" Message

**Cause:** The interactive graph visualization requires `pyvis`.

**Fix:**
```bash
pip install pyvis
```

The Knowledge Graph tab will fall back to a text-based display if pyvis is not available.

---

## 14. Quick Reference Card

Print this page and keep it at the podium during the demo.

### Services

| Service | URL | Purpose |
|---|---|---|
| Streamlit UI | `http://localhost:8521` | Main demo interface |
| FastAPI API | `http://localhost:8522` | REST API + Swagger docs |
| FastAPI Docs | `http://localhost:8522/docs` | Interactive API documentation |
| Milvus | `localhost:19530` | Vector database |

### Collections (11 total, 3,567,622 vectors)

| Collection | Vectors | Badge Color |
|---|---|---|
| cart_literature | 5,047 | Blue |
| cart_trials | 973 | Green |
| cart_constructs | 41 | Purple |
| cart_assays | 75 | Yellow |
| cart_manufacturing | 56 | Orange |
| cart_safety | 71 | Red |
| cart_biomarkers | 60 | Teal |
| cart_regulatory | 40 | Indigo |
| cart_sequences | 40 | Pink |
| cart_realworld | 54 | Brown |
| genomic_evidence | 3,561,170 | Teal |

### Demo Queries -- In Order

| # | Query | What It Shows | Time |
|---|---|---|---|
| 1 | Why do CD19 CAR-T therapies fail in relapsed B-ALL? | Basic RAG, citations, evidence panel | 5 min |
| 2 | Compare 4-1BB vs CD28 costimulatory domains | Comparative analysis, dual retrieval, structured tables | 5 min |
| 3 | What manufacturing parameters predict clinical response? | Cross-functional, manufacturing + clinical + biomarker | 5 min |
| 4 | Which biomarkers best predict CRS severity? | Safety + biomarker collections, toxicity knowledge | 3 min |
| 5 | What genomic variants in CD19 or BCMA pathway genes affect CAR-T response? | Genomic bridge, 3.5M variants, precision medicine | 3 min |
| 6 | What patents cover bispecific CAR-T constructs targeting CD19 and CD22? | Patent search, IP landscape | 2 min |
| 7 | How does scFv humanization reduce immunogenicity risk in CAR-T therapy? | Immunogenicity, HLA, sequence data | 2 min |

### All 13 Built-In Demo Query Buttons

These appear in the sidebar under "Demo Queries":

1. Why do CD19 CAR-T therapies fail in relapsed B-ALL?
2. Compare 4-1BB vs CD28 costimulatory domains
3. What manufacturing parameters predict response?
4. BCMA CAR-T resistance mechanisms in myeloma
5. How does T-cell exhaustion affect persistence?
6. What are the long-term safety signals for CD19 CAR-T products?
7. Which biomarkers best predict CRS severity?
8. Compare the FDA regulatory pathway of Kymriah vs Yescarta
9. What is the binding affinity of FMC63 scFv?
10. How do real-world CAR-T outcomes compare between academic and community centers?
11. What genomic variants in CD19 or BCMA pathway genes affect CAR-T response?
12. What patents cover bispecific CAR-T constructs targeting CD19 and CD22?
13. How does scFv humanization reduce immunogenicity risk in CAR-T therapy?

### API Endpoints Quick Reference

```bash
# Health check
curl http://localhost:8522/health

# Collection stats
curl http://localhost:8522/collections

# Evidence-only search (fast, no LLM)
curl -X POST http://localhost:8522/search \
  -H "Content-Type: application/json" \
  -d '{"question": "YOUR QUERY HERE"}'

# Full RAG query (retrieve + LLM)
curl -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{"question": "YOUR QUERY HERE"}'

# Cross-collection entity search
curl -X POST http://localhost:8522/find-related \
  -H "Content-Type: application/json" \
  -d '{"entity": "Yescarta", "top_k": 3}'

# Knowledge graph stats
curl http://localhost:8522/knowledge/stats

# Prometheus metrics
curl http://localhost:8522/metrics
```

### Key Numbers to Remember

| Metric | Value |
|---|---|
| Total vectors | **3,567,622** |
| Collections | **11** |
| Knowledge graph targets | **34** antigens |
| FDA-approved products | **6** (Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti) |
| Toxicity profiles | **17** |
| Manufacturing processes | **20** |
| Biomarkers | **23** |
| Immunogenicity topics | **6** |
| Query expansion keywords | **229** mapping to **1,961** terms |
| Retrieval latency | **20-30 ms** (11 collections in parallel) |
| Full RAG query time | **~24 seconds** (dominated by LLM) |
| Embedding model | BGE-small-en-v1.5 (384-dim) |
| LLM | Claude Sonnet 4.6 (Anthropic) |
| Hardware target | NVIDIA DGX Spark ($4,699) |
| License | Apache 2.0, open-source |

### Sidebar Controls Summary

| Control | Location | Purpose |
|---|---|---|
| Deep Research Mode | Top of sidebar | Toggles autonomous agent with sub-question decomposition |
| Target Antigen Filter | Below mode toggle | Filters results to a specific antigen (CD19, BCMA, etc.) |
| Development Stage | Below target filter | Filters by development stage (Target ID, CAR Design, etc.) |
| Date Range (From/To) | Middle of sidebar | Year-based filtering on publications and trials |
| Apply Date Filter | Below date range | Activates the date range filter |
| Collection Checkboxes | Lower sidebar | Enable/disable specific collections for search |
| Demo Queries | Bottom of sidebar | One-click buttons for 13 pre-built demo queries |

### Tabs Summary

| Tab | Purpose |
|---|---|
| **Chat** | Main query interface with streaming responses and evidence panel |
| **Knowledge Graph** | Interactive visualization of entity relationships (requires pyvis) |
| **Image Analysis** | Upload slides/images for claim extraction and evidence verification |

---

*Generated by HCLS AI Factory -- CAR-T Intelligence Agent v2.0.0 | Apache 2.0 | Adam Jones | March 2026*
