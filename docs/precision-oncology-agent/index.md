# Precision Oncology Agent

Closed-loop precision oncology clinical decision support -- from paired tumor-normal genomics to Molecular Tumor Board packets. Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory).

## Overview

The Precision Oncology Agent transforms raw genomic data (VCF files) into actionable clinical intelligence. It combines variant annotation, evidence retrieval, therapy ranking, trial matching, and outcomes learning into a closed-loop system that generates Molecular Tumor Board (MTB) packets for precision cancer treatment decisions.

| Collection | Records | Source |
|---|---|---|
| **Variants** | ~300 | CIViC, curated actionable variants |
| **Literature** | ~500 | PubMed oncology research |
| **Therapies** | ~120 | FDA-approved targeted/immuno/chemo |
| **Guidelines** | ~100 | NCCN, ASCO, ESMO guidelines |
| **Trials** | ~200 | ClinicalTrials.gov oncology trials |
| **Biomarkers** | ~80 | TMB, MSI-H, PD-L1, HRD, fusion panels |
| **Resistance** | ~80 | Resistance mechanisms and bypasses |
| **Pathways** | ~50 | Oncogenic signaling pathways |
| **Outcomes** | ~50 | Treatment outcomes (synthetic) |
| **Cases** | ~10 | Case snapshots (synthetic) |
| **Genomic Evidence** | *(read-only)* | Shared from Stage 2 RAG (3.5M vectors) |
| **Total** | **~1,490+ vectors** | **11 collections (10 owned + 1 shared)** |

### Example Queries

```
"What therapies target BRAF V600E in melanoma?"
"Compare EGFR TKI generations for NSCLC"
"Resistance mechanisms to osimertinib"
"NCCN recommendations for HER2+ breast cancer"
"Match clinical trials for MSI-H colorectal cancer"
"What is the role of TMB as a predictive biomarker for immunotherapy?"
```

All queries return grounded, cross-collection answers with clickable [PubMed](https://pubmed.ncbi.nlm.nih.gov/) and [ClinicalTrials.gov](https://clinicaltrials.gov/) citations.

### Comparative Analysis Mode

Comparative queries are **auto-detected** and produce structured side-by-side analysis with markdown tables, efficacy data, safety profiles, and clinical context.

```
"Compare osimertinib vs erlotinib for EGFR-mutant NSCLC"     -> TKI generation comparison
"BRAF+MEK inhibition vs immunotherapy for melanoma"           -> Modality comparison
"Pembrolizumab versus nivolumab"                              -> Product vs product
"Compare PARP inhibitors for BRCA-mutant ovarian cancer"      -> Drug class comparison
```

**How it works:** The engine detects "vs/versus/compare/difference between" keywords, parses two entities, resolves each against the knowledge graph (~40 actionable targets, ~30 therapy mappings, ~20 resistance mechanisms), runs **dual retrievals** with per-entity filtering, identifies shared/head-to-head evidence, and builds a comparative prompt instructing Claude to produce structured tables with efficacy, safety, biomarker, and guideline data.

| Feature | Detail |
|---|---|
| Entity types | Genes, drugs, drug classes, biomarkers, cancer types, pathways |
| Entity resolution | ~40 gene targets + ~30 therapy mappings + 50+ aliases |
| Dual retrieval | ~400 ms for 46 results (24 + 22 per entity) |
| Shared evidence | Head-to-head trials identified and highlighted separately |
| Structured output | Comparison table, MoA differences, efficacy, safety, guideline recs |
| Fallback | Unrecognized entities gracefully fall back to normal query path |

### MTB Packet Generation

The Precision Oncology Agent generates structured Molecular Tumor Board packets from patient data:

1. **VCF Upload** -- raw VCF text is parsed, extracting PASS variants with gene, consequence, and position
2. **Variant Annotation** -- each variant is classified against ~40 actionable targets using AMP/ASCO/CAP evidence tiers (A-D)
3. **Evidence Lookup** -- RAG retrieval for each actionable variant across literature, therapies, and guidelines
4. **Therapy Ranking** -- evidence-level-sorted therapy recommendations with resistance flags and contraindication checks
5. **Trial Matching** -- hybrid deterministic + semantic search against oncology clinical trials
6. **Open Questions** -- VUS variants, missing biomarkers, and evidence gaps flagged for MTB discussion

The resulting MTB packet is exported as Markdown, JSON, PDF, or FHIR R4 DiagnosticReport Bundle.

## Architecture

```
VCF / Patient Data
    |
    v
[Case Manager] ---- VCF parsing, variant extraction, actionability classification
    |
    v
[Knowledge Graph Lookup]
(~40 actionable targets, ~30 therapies, ~20 resistance, ~10 pathways, ~15 biomarkers)
    |
    v
[Parallel 11-Collection RAG Search] --- BGE-small-en-v1.5 (384-dim)
    |               |              |           |           |
    v               v              v           v           v
 Variants      Literature     Therapies   Guidelines    Trials
  ~300           ~500           ~120         ~100        ~200
    |               |              |           |           |
 Biomarkers    Resistance     Pathways    Outcomes      Cases
   ~80            ~80           ~50          ~50          ~10
    |                                                      |
    +-------------- genomic_evidence (3.5M, read-only) ----+
    |
    v
[Query Expansion] (12 maps, ~120 keywords -> ~700 terms)
    |
    v
[Evidence Synthesis]
    |
    +--- [Therapy Ranker] -- evidence-level sort, resistance check, contraindication
    +--- [Trial Matcher] --- deterministic filter + semantic search + composite scoring
    +--- [Cross-Modal] ----- variant severity -> imaging, variant actionability -> drug discovery
    |
    v
[Claude Sonnet 4.6 LLM] -> Grounded response with citations
    |
    v
[MTB Packet / Clinical Report]
    |
    +--- Markdown report
    +--- JSON export
    +--- PDF (NVIDIA-themed, ReportLab)
    +--- FHIR R4 DiagnosticReport Bundle (SNOMED CT, LOINC coded)
```

Built on the HCLS AI Factory platform:

- **Vector DB:** Milvus 2.4 with IVF_FLAT/COSINE indexes (nlist=1024, nprobe=16)
- **Embeddings:** BGE-small-en-v1.5 (384-dim)
- **LLM:** Claude Sonnet 4.6 (Anthropic API)
- **UI:** Streamlit MTB Workbench (port 8526)
- **API:** FastAPI REST server (port 8527)
- **Hardware target:** NVIDIA DGX Spark ($3,999)

## Setup

### Prerequisites

- Python 3.10+
- Milvus 2.4 running on `localhost:19530`
- `ANTHROPIC_API_KEY` environment variable (or in `rag-chat-pipeline/.env`)

### Install

```bash
cd ai_agent_adds/precision_oncology_agent
pip install -r requirements.txt
```

### 1. Create Collections and Seed Data

```bash
python3 scripts/setup_collections.py --seed
```

This creates 11 Milvus collections (10 owned + 1 read-only) with IVF_FLAT indexes and seeds actionable variants, therapies, guidelines, resistance mechanisms, pathways, and biomarker panels.

### 2. Ingest PubMed Literature (~15 min)

```bash
python3 scripts/ingest_pubmed.py --max-results 5000
```

Fetches oncology abstracts via NCBI E-utilities, classifies by cancer type and gene, embeds with BGE-small, and stores in `onco_literature`.

### 3. Ingest Clinical Trials (~3 min)

```bash
python3 scripts/ingest_clinical_trials.py --max-results 1500
```

Fetches oncology trials via ClinicalTrials.gov API v2, extracts biomarker criteria, embeds, and stores in `onco_trials`.

### 4. Ingest CIViC Variants (~2 min)

```bash
python3 scripts/ingest_civic.py
```

Fetches clinically actionable variants from the CIViC database, maps evidence levels to AMP/ASCO/CAP tiers, and stores in `onco_variants`.

### 5. Validate

```bash
python3 scripts/validate_e2e.py
```

Runs end-to-end tests: collection stats, single-collection search, multi-collection `search_all()`, filtered search, and demo queries.

### 6. Launch UI

```bash
streamlit run app/oncology_ui.py --server.port 8526
```

### 7. Launch API

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8527
```

## Project Structure

```
precision_oncology_agent/agent/
├── src/
│   ├── models.py                  # Pydantic data models (497 lines)
│   ├── collections.py             # 11 Milvus collection schemas + manager (606 lines)
│   ├── knowledge.py               # Knowledge graph: targets, therapies, resistance (1,194 lines)
│   ├── query_expansion.py         # 12 expansion maps (676 lines)
│   ├── rag_engine.py              # Multi-collection RAG + comparative analysis (780 lines)
│   ├── agent.py                   # Plan-search-synthesize pipeline (489 lines)
│   ├── case_manager.py            # VCF parsing + MTB packet generation (509 lines)
│   ├── trial_matcher.py           # Hybrid deterministic + semantic matching (393 lines)
│   ├── therapy_ranker.py          # Evidence-based therapy ranking (552 lines)
│   ├── cross_modal.py             # Cross-modal triggers to imaging + drug discovery (395 lines)
│   ├── export.py                  # Markdown, JSON, PDF, FHIR R4 export (876 lines)
│   ├── metrics.py                 # Prometheus metrics (362 lines)
│   ├── scheduler.py               # Data ingestion scheduler (263 lines)
│   ├── ingest/
│   │   ├── base.py                # Base ingest pipeline (249 lines)
│   │   ├── civic_parser.py        # CIViC actionable variant ingest (340 lines)
│   │   ├── oncokb_parser.py       # OncoKB data parser (104 lines)
│   │   ├── literature_parser.py   # PubMed E-utilities ingest (248 lines)
│   │   ├── clinical_trials_parser.py  # ClinicalTrials.gov API v2 (279 lines)
│   │   ├── guideline_parser.py    # NCCN/ASCO/ESMO guideline parser (168 lines)
│   │   ├── pathway_parser.py      # Signaling pathway parser (121 lines)
│   │   ├── resistance_parser.py   # Resistance mechanism parser (125 lines)
│   │   └── outcome_parser.py      # Outcome record parser (158 lines)
│   └── utils/
│       ├── vcf_parser.py          # VCF file parsing utilities (361 lines)
│       └── pubmed_client.py       # NCBI E-utilities HTTP client (296 lines)
├── app/
│   └── oncology_ui.py             # Streamlit MTB Workbench (703 lines)
├── api/
│   ├── main.py                    # FastAPI REST server (347 lines)
│   └── routes/
│       ├── meta_agent.py          # /api/ask, /api/deep-research (169 lines)
│       ├── cases.py               # /api/cases, /api/cases/{id}/mtb (234 lines)
│       ├── trials.py              # /api/trials/match (153 lines)
│       ├── reports.py             # /api/reports/{format} (236 lines)
│       └── events.py              # /api/events, cross-modal triggers (89 lines)
├── config/
│   └── settings.py                # Pydantic BaseSettings (109 lines)
├── tests/
│   └── conftest.py                # Test fixtures (214 lines)
├── requirements.txt
└── LICENSE                        # Apache 2.0
```

**63 Python files | ~17,855 lines of code | Apache 2.0**

## Knowledge Graph

| Category | Count | Examples |
|---|---|---|
| Actionable Targets | ~40 | BRAF, EGFR, ALK, ROS1, KRAS G12C, HER2, NTRK, RET, MET, FGFR, PIK3CA, BRCA1/2, IDH1/2, ESR1, TP53, PTEN... |
| Therapy Mappings | ~30 | vemurafenib, osimertinib, pembrolizumab, sotorasib, lorlatinib, olaparib, trastuzumab deruxtecan... |
| Resistance Mechanisms | ~20 | EGFR T790M, EGFR C797S, MET amplification bypass, BRAF amplification, NRAS activation, BRCA reversion... |
| Signaling Pathways | ~10 | MAPK, PI3K/AKT/mTOR, DDR, cell cycle, apoptosis, Wnt, Notch, Hedgehog, JAK/STAT, angiogenesis |
| Biomarker Panels | ~15 | TMB-H, MSI-H/dMMR, PD-L1 TPS/CPS, HRD, NTRK fusion, ALK rearrangement, ROS1 fusion, ctDNA, FGFR |
| Entity Aliases | ~50+ | Cancer type aliases (lung->NSCLC, CRC->COLORECTAL), drug brand names, gene synonyms |

## Performance

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores):

| Metric | Value |
|---|---|
| Single-collection search | < 200 ms |
| Cross-collection RAG query (11 collections) | < 5 s |
| MTB packet generation (full workflow) | < 30 s |
| Trial matching (deterministic + semantic) | < 10 s |
| Therapy ranking (variant + biomarker driven) | < 5 s |
| Comparative dual retrieval | ~400 ms |
| Full RAG query (search + Claude) | ~24 s |
| Cosine similarity scores | 0.72 - 0.92 |
| Embedding dimension | 384 (BGE-small-en-v1.5) |

## Service Ports

| Service | Port | Protocol |
|---|---|---|
| Streamlit MTB Workbench | 8526 | HTTP |
| FastAPI REST Server | 8527 | HTTP |
| Milvus gRPC | 19530 | gRPC |

## Credits

- **Adam Jones**
- **Apache 2.0 License**
