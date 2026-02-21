# CAR-T Intelligence Agent

Cross-functional intelligence across the CAR-T cell therapy development lifecycle. Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory).

## Overview

The CAR-T Intelligence Agent breaks down data silos across the 5 stages of CAR-T development. It searches across **all data sources simultaneously** and synthesizes cross-functional insights powered by Claude.

| Collection | Records | Source |
|---|---|---|
| **Literature** | 5,047 | PubMed abstracts via NCBI E-utilities |
| **Clinical Trials** | 973 | ClinicalTrials.gov API v2 |
| **CAR Constructs** | 6 | 6 FDA-approved CAR-T products |
| **Assay Results** | 45 | Curated from landmark papers (ELIANA, ZUMA-1, KarMMa, CARTITUDE-1, etc.) |
| **Manufacturing** | 30 | Curated CMC/process data (transduction, expansion, release, cryo, logistics) |
| **Safety** | — | Pharmacovigilance, CRS/ICANS profiles |
| **Biomarkers** | — | CRS prediction, exhaustion monitoring |
| **Regulatory** | — | Approval timelines, post-marketing requirements |
| **Sequences** | — | Molecular binding, scFv sequences |
| **Real-World Evidence** | — | Registry outcomes, real-world data |
| **Genomic Evidence** | *(read-only)* | Shared from Stage 2 RAG pipeline (Milvus) |
| **Total** | **6,266+ vectors** | **11 collections (10 owned + 1 read-only)** |

### Example Queries

```
"Why do CD19 CAR-T therapies fail in relapsed B-ALL?"
"Compare 4-1BB vs CD28 costimulatory domains for DLBCL"
"What manufacturing parameters predict clinical response?"
"BCMA CAR-T resistance mechanisms in multiple myeloma"
"How does T-cell exhaustion affect CAR-T persistence?"
```

All queries return grounded, cross-collection answers with clickable [Literature:PMID](https://pubmed.ncbi.nlm.nih.gov/) and [Trial:NCT...](https://clinicaltrials.gov/) citations.

### Comparative Analysis Mode

Comparative queries are **auto-detected** and produce structured side-by-side analysis with markdown tables, advantages/limitations, and clinical context.

```
"Compare CD19 vs BCMA"                              → Target vs target
"Compare 4-1BB vs CD28 costimulatory domains"        → Costimulatory domain comparison
"Kymriah versus Carvykti"                            → Product vs product (resolves to CD19 vs BCMA)
"Compare CRS and ICANS toxicity"                     → Toxicity profile comparison
```

**How it works:** The engine detects "vs/versus/compare" keywords, parses two entities, resolves each against the knowledge graph (25 antigens, 6 products, 8 toxicities, 10 manufacturing processes), runs **dual retrievals** with per-entity filtering, and builds a comparative prompt that instructs Claude to produce structured tables. The evidence panel groups results by entity with color-coded headers.

| Feature | Detail |
|---|---|
| Entity types | Targets, FDA products, costimulatory domains, toxicities, manufacturing |
| Entity resolution | 25 antigens + 39+ product/domain/biomarker/regulatory aliases |
| Dual retrieval | ~365ms for 46 results (24 + 22 per entity) |
| Structured output | Comparison table, advantages, limitations, clinical context |
| Fallback | Unrecognized entities gracefully fall back to normal query path |

## Architecture

```
User Query
    |
    v
[Comparative Detection] ──── "X vs Y" detected? ──── YES ──┐
    |                                                        |
    NO                                              [Parse Two Entities]
    |                                              (resolve via knowledge graph)
    v                                                        |
[BGE-small-en-v1.5 Embedding]                      [Dual Retrieval]
(384-dim, asymmetric query prefix)                  (Entity A + Entity B)
    |                                                        |
    v                                                        v
[Parallel Search: 11 Milvus Collections]     [Comparative Prompt Builder]
(IVF_FLAT / COSINE)                         (tables + pros/cons format)
    |               |           |                            |
    v               v           v                            |
Literature      Trials     Constructs                        |
 5,047           973           6                             |
    |               |           |                            |
 Assays      Manufacturing     |                             |
   45             30           |                             |
    |               |           |                            |
    +-------+-------+----------+                             |
            |                                                |
            v                                                v
    [Query Expansion] (12 maps, 169 keywords -> 1,496 terms)
            |
            v
    [Knowledge Graph Augmentation]
    (25 antigens, 8 toxicities, 10 mfg processes)
            |
            v
    [Claude LLM] -> Grounded response with citations
```

Built on the HCLS AI Factory platform:

- **Vector DB:** Milvus 2.4 with IVF_FLAT/COSINE indexes (nlist=1024, nprobe=16)
- **Embeddings:** BGE-small-en-v1.5 (384-dim)
- **LLM:** Claude Sonnet 4.6 (Anthropic API)
- **UI:** Streamlit (port 8521)
- **Hardware target:** NVIDIA DGX Spark ($3,999)

## Setup

### Prerequisites

- Python 3.10+
- Milvus 2.4 running on `localhost:19530`
- `ANTHROPIC_API_KEY` environment variable (or in `rag-chat-pipeline/.env`)

### Install

```bash
cd ai_agent_adds/cart_intelligence_agent
pip install -r requirements.txt
```

### 1. Create Collections and Seed FDA Constructs

```bash
python3 scripts/setup_collections.py --seed-constructs
```

This creates 11 Milvus collections (10 owned + 1 read-only) with IVF_FLAT indexes and inserts 6 FDA-approved CAR-T products (Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti).

### 2. Ingest PubMed Literature (~15 min)

```bash
python3 scripts/ingest_pubmed.py --max-results 5000
```

Fetches CAR-T abstracts via NCBI E-utilities (esearch + efetch), classifies by development stage, extracts target antigens, embeds with BGE-small, and stores in `cart_literature`.

### 3. Ingest Clinical Trials (~3 min)

```bash
python3 scripts/ingest_clinical_trials.py --max-results 1500
```

Fetches CAR-T trials via ClinicalTrials.gov API v2, extracts phase/status/sponsor/antigen/generation, embeds, and stores in `cart_trials`.

### 4. Seed Assay Data (~30 sec)

```bash
python3 scripts/seed_assays.py
```

Inserts 45 curated assay records from landmark CAR-T publications (ELIANA, ZUMA-1, KarMMa, CARTITUDE-1, etc.) covering cytotoxicity, cytokine, persistence, exhaustion, and resistance data.

### 5. Seed Manufacturing Data (~30 sec)

```bash
python3 scripts/seed_manufacturing.py
```

Inserts 30 curated manufacturing/CMC records covering transduction, expansion, harvest, cryopreservation, release testing, logistics, cost analysis, and emerging platforms (POC, allogeneic, non-viral).

### 6. Validate

```bash
python3 scripts/validate_e2e.py
```

Runs 5 tests: collection stats, single-collection search, multi-collection `search_all()`, filtered search (`target_antigen == "CD19"`), and all demo queries.

### 7. Run Integration Test (requires API key)

```bash
python3 scripts/test_rag_pipeline.py
```

Tests the full RAG pipeline: embed -> search_all -> knowledge graph -> Claude LLM response generation. Validates both synchronous and streaming modes.

### 8. Launch UI

```bash
streamlit run app/cart_ui.py --server.port 8521
```

## Project Structure

```
cart_intelligence_agent/
├── Docs/
│   └── CART_Intelligence_Agent_Design.md  # Architecture design document
├── src/
│   ├── models.py                  # Pydantic data models (16 models + enums)
│   ├── collections.py             # 11 Milvus collection schemas + manager
│   ├── knowledge.py               # Knowledge graph (25 targets, 8 toxicities, 10 mfg)
│   ├── query_expansion.py         # 12 expansion maps (169 keywords -> 1,496 terms)
│   ├── rag_engine.py              # Multi-collection RAG engine + comparative analysis + Claude
│   ├── agent.py                   # CAR-T Intelligence Agent (plan -> search -> synthesize)
│   ├── ingest/
│   │   ├── base.py                # Base ingest pipeline (fetch -> parse -> embed -> store)
│   │   ├── literature_parser.py   # PubMed NCBI E-utilities ingest
│   │   ├── clinical_trials_parser.py  # ClinicalTrials.gov API v2 ingest
│   │   ├── construct_parser.py    # CAR construct data parser
│   │   ├── assay_parser.py        # Assay data parser
│   │   └── manufacturing_parser.py # Manufacturing/CMC data parser
│   └── utils/
│       └── pubmed_client.py       # NCBI E-utilities HTTP client
├── app/
│   └── cart_ui.py                 # Streamlit chat interface (NVIDIA theme, comparative mode)
├── config/
│   └── settings.py                # Pydantic BaseSettings configuration
├── data/
│   └── reference/
│       ├── assay_seed_data.json   # 45 curated assay records from landmark papers
│       └── manufacturing_seed_data.json  # 30 curated manufacturing/CMC records
├── scripts/
│   ├── setup_collections.py       # Create collections + seed FDA constructs
│   ├── ingest_pubmed.py           # CLI: ingest PubMed CAR-T literature
│   ├── ingest_clinical_trials.py  # CLI: ingest ClinicalTrials.gov trials
│   ├── seed_assays.py             # Seed assay data from published papers
│   ├── seed_manufacturing.py      # Seed manufacturing/CMC data
│   ├── validate_e2e.py            # End-to-end data layer validation
│   ├── test_rag_pipeline.py       # Full RAG + LLM integration test
│   └── seed_knowledge.py          # Export knowledge graph to JSON
├── requirements.txt
└── LICENSE                        # Apache 2.0
```

**55 Python files | ~16,748 lines | Apache 2.0**

## Knowledge Graph

| Component | Count | Examples |
|---|---|---|
| Target Antigens | 25 | CD19, BCMA, CD22, CD20, CD30, HER2, GPC3, EGFR, Mesothelin, GPRC5D, ... |
| FDA-Approved Products | 6 | Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti |
| Entity Aliases | 39+ | Product names, generic names, costimulatory domains (for comparative resolution) |
| Toxicity Profiles | 8 | CRS, ICANS, B-cell aplasia, HLH/MAS, cytopenias, TLS, GvHD, on-target/off-tumor |
| Manufacturing Processes | 10 | Transduction, expansion, leukapheresis, cryopreservation, release testing, ... |
| Biomarkers | 15 | CRS prediction, T-cell exhaustion, persistence, cytokine panels, ... |
| Regulatory Histories | 6 | Approval timelines, post-marketing requirements for all FDA products |
| Immunogenicity Topics | 6 | ADA, immunogenicity assays, risk factors, ... |
| Query Expansion Maps | 12 | Target Antigen, Disease, Toxicity, Manufacturing, Mechanism, Construct, Safety, Biomarker, Regulatory, Sequence, RealWorld, Immunogenicity |
| Expansion Keywords | 169 | Mapping to 1,496 related terms |

## Performance

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified memory):

| Metric | Value |
|---|---|
| PubMed ingest (5,047 abstracts) | ~15 min |
| ClinicalTrials.gov ingest (973 trials) | ~3 min |
| Assay seed ingest (45 records) | ~30 sec |
| Manufacturing seed ingest (30 records) | ~30 sec |
| Vector search (11 collections, top-5 each) | 12-16 ms (cached) |
| Comparative dual retrieval (2x11 collections) | ~365 ms |
| Full RAG query (search + Claude) | ~24 sec |
| Comparative RAG query (dual search + Claude) | ~30 sec |
| Cosine similarity scores | 0.74 - 0.90 |

## Status

- **Week 1 (Scaffold)** -- Complete. Architecture, data models, collection schemas, knowledge graph, ingest pipelines, RAG engine, agent, and Streamlit UI.
- **Week 2 Days 1-3 (Data)** -- Complete. PubMed (5,047) + ClinicalTrials.gov (973) + FDA constructs (6) ingested. End-to-end validation passing.
- **Week 2 Days 4-5 (Integration)** -- Complete. Full RAG pipeline with Claude LLM generating grounded cross-functional answers. Streamlit UI working.
- **Week 2 Day 5+ (Assay + Manufacturing Data)** -- Complete. 45 curated assay records + 30 manufacturing/CMC records seeded. All 11 collections populated (10 owned + 1 read-only). Total: 6,266 owned vectors.
- **Week 3 (UI + Analysis)** -- Complete. Clickable PubMed/ClinicalTrials.gov citation links, collapsible evidence panel with collection badges, and **Comparative Analysis Mode** with auto-detection, dual retrieval, entity-grouped evidence, and structured markdown tables.

## Credits

- **Adam Jones** -- HCLS AI Factory, 14+ years genomic research
- **Apache 2.0 License**
