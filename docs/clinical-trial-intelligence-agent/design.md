# Clinical Trial Intelligence Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## 1. Executive Summary

The Clinical Trial Intelligence Agent extends the HCLS AI Factory platform to deliver RAG-powered decision support across the entire clinical trial lifecycle. It integrates protocol design optimization, patient-trial matching, site selection, safety signal detection, and competitive landscape analysis into a single intelligence platform, enabling pharmaceutical R&D teams to make evidence-based decisions faster and with greater confidence.

The platform enables cross-functional queries like *"Design an adaptive Phase II/III protocol for EGFR-mutant NSCLC with a biomarker-enriched population"* that simultaneously search protocol templates, eligibility criteria, endpoint definitions, regulatory precedents, safety databases, and competitive intelligence — returning grounded recommendations with citations to landmark trials and regulatory guidance.

### Key Results

| Metric | Value |
|---|---|
| Total vectors indexed | **~251,500** across 14 Milvus collections (13 owned + 1 read-only) |
| Clinical workflows | **10** (protocol design, patient matching, site selection, eligibility, adaptive, safety signal, regulatory, competitive, diversity, DCT) |
| Decision support engines | **5** + 1 Historical Success Estimator |
| Landmark trials in knowledge base | **40** across 13 therapeutic areas |
| Regulatory agencies modeled | **9** (FDA, EMA, PMDA, NMPA, Health Canada, TGA, MHRA, Swissmedic, ANVISA) |
| Entity aliases | **140** query expansion mappings |
| API endpoints | **26** |
| Lines of code | **22,607** |
| Test suite | **769** tests (100% pass, 0.47s runtime) |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Clinical Trial Agent Role |
|---|---|
| **DataStore** | Raw files: ClinicalTrials.gov JSON, PubMed XML, FDA regulatory docs, site performance data |
| **DataEngine** | Event-driven ingest pipelines with 7+ parsers (protocol, eligibility, endpoint, site, safety, regulatory, literature) |
| **DataBase** | 14 Milvus collections (13 owned + 1 read-only) + knowledge base (40 trials, 13 areas, 9 agencies) |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + query expansion (10 maps, 140 aliases) |
| **AgentEngine** | ClinicalTrialAgent (Plan-Search-Evaluate-Synthesize) + Streamlit UI (5 tabs) |

### 2.2 System Diagram

```
                        ┌─────────────────────────────────┐
                        │    Streamlit Chat UI (8128)       │
                        │    5 tabs: Intelligence |         │
                        │    Matching | Protocol |          │
                        │    Competitive | Dashboard        │
                        └──────────────┬──────────────────┘
                                       │
                        ┌──────────────▼──────────────────┐
                        │     FastAPI REST API (8538)       │
                        │     26 endpoints, CORS, Auth      │
                        │     Rate limiting, Metrics         │
                        └──────────────┬──────────────────┘
                                       │
                        ┌──────────────▼──────────────────┐
                        │     ClinicalTrialAgent            │
                        │  Plan → Search → Evaluate →       │
                        │  Synthesize                        │
                        └──────────────┬──────────────────┘
                                       │
        ┌──────────────────────────────┼──────────────────────────────┐
        │                              │                              │
┌───────▼────────┐          ┌──────────▼──────────┐       ┌──────────▼──────────┐
│ Workflows (10) │          │ Decision Engines (5) │       │ RAG Engine           │
│                │          │                      │       │                      │
│ Protocol Design│          │ Confidence Scorer    │       │ BGE-small-en-v1.5    │
│ Patient Match  │          │ Complexity Estimator │       │ (384-dim embedding)  │
│ Site Selection │          │ Enrollment Predictor │       │         │            │
│ Eligibility    │          │ Eligibility Scorer   │       │         ▼            │
│ Adaptive Design│          │ Competitive Ranker   │       │ Parallel Search      │
│ Safety Signal  │          │                      │       │ 14 Milvus Collections│
│ Regulatory     │          │ + Historical Success │       │ (ThreadPoolExecutor) │
│ Competitive    │          │   Estimator          │       │         │            │
│ Diversity      │          │                      │       │         ▼            │
│ DCT Planning   │          │                      │       │ Claude Sonnet 4.6    │
└───────┬────────┘          └──────────┬──────────┘       └──────────────────────┘
        │                              │
┌───────▼──────────────────────────────▼──────────────────────────────┐
│                  Milvus 2.4 — 14 Collections                        │
│                                                                      │
│  trial_protocols (5K)       trial_eligibility (50K)                  │
│  trial_endpoints (20K)      trial_sites (30K)                        │
│  trial_investigators (5K)   trial_results (3K)                       │
│  trial_regulatory (2K)      trial_literature (10K)                   │
│  trial_biomarkers (3K)      trial_safety (20K)                       │
│  trial_rwe (2K)             trial_adaptive (500)                     │
│  trial_guidelines (1K)      genomic_evidence (100K) [shared]         │
└──────────────────────────────────────────────────────────────────────┘
```

---

## 3. Data Collections — Actual State

All 14 collections (13 owned + 1 shared read-only) are populated and searchable.

### 3.1 Collection Catalog

| # | Collection | Est. Records | Weight | Primary Use |
|---|---|---|---|---|
| 1 | `trial_protocols` | 5,000 | 0.10 | Protocol design, competitive intelligence |
| 2 | `trial_eligibility` | 50,000 | 0.09 | Patient matching, eligibility optimization |
| 3 | `trial_endpoints` | 20,000 | 0.08 | Protocol design, adaptive design evaluation |
| 4 | `trial_sites` | 30,000 | 0.07 | Site selection, diversity assessment |
| 5 | `trial_investigators` | 5,000 | 0.05 | Site selection, competitive intelligence |
| 6 | `trial_results` | 3,000 | 0.09 | Protocol design, competitive intelligence |
| 7 | `trial_regulatory` | 2,000 | 0.07 | Regulatory document generation |
| 8 | `trial_literature` | 10,000 | 0.08 | Evidence synthesis, protocol design |
| 9 | `trial_biomarkers` | 3,000 | 0.07 | Patient matching, biomarker strategy |
| 10 | `trial_safety` | 20,000 | 0.08 | Safety signal detection, regulatory docs |
| 11 | `trial_rwe` | 2,000 | 0.06 | Eligibility optimization, diversity planning |
| 12 | `trial_adaptive` | 500 | 0.05 | Adaptive design evaluation |
| 13 | `trial_guidelines` | 1,000 | 0.08 | All workflows (regulatory reference) |
| 14 | `genomic_evidence` | 100,000 | 0.03 | Cross-modal genomic queries |

### 3.2 Index Configuration (all collections)

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 1024 (protocols, eligibility, safety), 256 (others) |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Knowledge Base

### 4.1 Landmark Trials (40 entries)

Each entry includes: NCT ID, trial name, therapeutic area, phase, design type, primary endpoint, key result, regulatory outcome, and lessons learned.

| Therapeutic Area | Count | Example Trials |
|---|---|---|
| **Oncology** | 12 | KEYNOTE-024, CheckMate-067, DESTINY-Breast03, ELIANA |
| **Cardiology** | 5 | DAPA-HF, EMPEROR-Reduced, PARADIGM-HF |
| **Neurology** | 4 | CLARITY AD, TRAILBLAZER-ALZ, EMERGE/ENGAGE |
| **Immunology** | 4 | TRANSFORM, RINVOQ, SKYRIZI |
| **Rare Disease** | 3 | FIREFISH, SUNFISH, SPRINT |
| **Infectious Disease** | 3 | RECOVERY, SOLIDARITY, COVE |
| **Metabolic** | 3 | SURPASS, SELECT, STEP |
| **Other** | 6 | Various therapeutic areas |

### 4.2 Regulatory Agencies (9 entries)

| Agency | Jurisdiction | Key Pathways |
|---|---|---|
| **FDA** | United States | Breakthrough Therapy, Accelerated Approval, RMAT, Fast Track |
| **EMA** | European Union | PRIME, Conditional MA, Orphan Designation |
| **PMDA** | Japan | SAKIGAKE, Conditional/Time-Limited Approval |
| **NMPA** | China | Priority Review, Breakthrough Therapy |
| **Health Canada** | Canada | Priority Review, NOC/c |
| **TGA** | Australia | Priority Review, Provisional Approval |
| **MHRA** | United Kingdom | ILAP, Conditional MA |
| **Swissmedic** | Switzerland | Fast-track, Temporary Authorization |
| **ANVISA** | Brazil | Priority Review |

### 4.3 Adaptive Design Templates (9 entries)

| Design Type | Use Case | Regulatory Precedent |
|---|---|---|
| Bayesian adaptive randomization | Biomarker-driven oncology | I-SPY 2 |
| Platform trial | Multi-arm, multi-stage | RECOVERY |
| Seamless Phase II/III | Dose selection + confirmatory | KEYNOTE-024 |
| Group sequential | Early stopping for efficacy/futility | Most pivotal trials |
| Response-adaptive | Rare disease, small populations | SUNFISH |
| Biomarker-adaptive | Enrichment based on interim | CheckMate-227 |
| Sample size re-estimation | Adaptive enrollment | PARADIGM-HF |
| Master protocol (basket) | Histology-independent | NCI-MATCH |
| Master protocol (umbrella) | Biomarker-stratified | Lung-MAP |

---

## 5. Clinical Workflows

### 5.1 Workflow Catalog

| # | Workflow | Clinical Question | Key Collections |
|---|---|---|---|
| 1 | Protocol Design | "What endpoints and design would optimize this Phase II/III?" | protocols, endpoints, results, adaptive |
| 2 | Patient Matching | "Which patients meet eligibility for this trial?" | eligibility, biomarkers, rwe |
| 3 | Site Selection | "Which sites will enroll fastest with the best data quality?" | sites, investigators, results |
| 4 | Eligibility Optimization | "How can we broaden eligibility without compromising safety?" | eligibility, safety, rwe, guidelines |
| 5 | Adaptive Design | "Should we use Bayesian adaptive randomization or group sequential?" | adaptive, endpoints, results |
| 6 | Safety Signal | "Are there emerging safety signals in this ongoing trial?" | safety, literature, regulatory |
| 7 | Regulatory Strategy | "What regulatory pathway should we pursue for accelerated approval?" | regulatory, guidelines, results |
| 8 | Competitive Landscape | "Who else is developing therapies for this indication?" | protocols, results, literature |
| 9 | Diversity Planning | "How do we ensure representative enrollment across demographics?" | sites, rwe, eligibility |
| 10 | DCT Planning | "Which trial activities can be decentralized?" | sites, guidelines, rwe |

### 5.2 Decision Support Engines

| Engine | Function | Output |
|---|---|---|
| **Confidence Scorer** | Assesses evidence strength for recommendations | 0-100 confidence score with supporting evidence count |
| **Complexity Estimator** | Evaluates protocol complexity against benchmarks | Complexity index, simplified alternatives |
| **Enrollment Predictor** | Projects enrollment timelines from site data | Months to full enrollment, at-risk sites |
| **Eligibility Scorer** | Rates eligibility criteria restrictiveness | Restrictiveness score, broadening suggestions |
| **Competitive Ranker** | Ranks competitive trials by threat level | Ranked competitor list with differentiation gaps |
| **Historical Success Estimator** | Predicts trial success probability from historical data | Success probability by phase, indication, design |

---

## 6. Multi-Collection RAG Engine

### 6.1 Search Flow

```
User Query: "Design a Phase II protocol for KRAS G12C NSCLC"
    │
    ├── 1. Embed query with BGE asymmetric prefix               [< 5 ms]
    │
    ├── 2. Parallel search across 14 collections (top-5 each)   [12-18 ms]
    │   ├── trial_protocols:    KRAS G12C NSCLC trials           (score: 0.82-0.90)
    │   ├── trial_results:      Sotorasib/adagrasib outcomes     (score: 0.78-0.88)
    │   ├── trial_endpoints:    NSCLC Phase II endpoints         (score: 0.75-0.85)
    │   ├── trial_eligibility:  KRAS-selected criteria           (score: 0.74-0.82)
    │   └── trial_adaptive:     Biomarker-adaptive designs       (score: 0.70-0.80)
    │
    ├── 3. Query expansion: "KRAS G12C NSCLC" →                 [< 1 ms]
    │      [KRAS, sotorasib, adagrasib, non-small cell,
    │       biomarker-selected, targeted therapy, ...]
    │
    ├── 4. Weighted merge + deduplicate (cap at 30 results)      [< 1 ms]
    │
    ├── 5. Knowledge base augmentation                           [< 1 ms]
    │      NSCLC → landmark trials, regulatory pathways,
    │      historical success rates, competitive landscape
    │
    └── 6. Stream Claude Sonnet 4.6 response                    [~22-26 sec]
           Grounded protocol recommendation with
           design rationale, endpoint selection,
           and regulatory pathway guidance
```

**Total: ~26 sec** (dominated by LLM generation; retrieval is ~25 ms)

### 6.2 Embedding Strategy

**Model:** BGE-small-en-v1.5 (BAAI)

| Mode | Prefix | Usage |
|---|---|---|
| **Query** | `"Represent this sentence for searching relevant passages: "` | User questions |
| **Document** | None (raw text) | Ingested records |

---

## 7. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 7.1 Search Performance

| Operation | Latency | Notes |
|---|---|---|
| Single collection search (top-5) | 3-5 ms | Milvus IVF_FLAT with cached index |
| 14-collection parallel search (top-5 each) | 12-18 ms | ThreadPoolExecutor, 70 total results |
| Query expansion + filtered search | 8-12 ms | Up to 5 expanded terms, applicable collections |
| Knowledge base augmentation | < 1 ms | In-memory dictionary lookup |
| Full retrieve() pipeline | 22-32 ms | Embed + search + expand + merge + knowledge |

### 7.2 RAG Query Performance

| Operation | Latency | Notes |
|---|---|---|
| Full query (retrieve + Claude generate) | ~26 sec | Dominated by LLM generation |
| Streaming query (time to first token) | ~3 sec | Evidence returned immediately |
| Response length | 1000-2500 chars | Grounded answer with citations |
| Token count | 500-1200 tokens | Claude Sonnet 4.6 output |

### 7.3 Decision Engine Performance

| Engine | Latency |
|---|---|
| Confidence Scorer | <20 ms |
| Complexity Estimator | <30 ms |
| Enrollment Predictor | <50 ms |
| Eligibility Scorer | <20 ms |
| Competitive Ranker | <40 ms |
| Historical Success Estimator | <30 ms |
| **All engines combined** | **<200 ms** |

---

## 8. Infrastructure

### 8.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4, localhost:19530 |
| Embeddings | BGE-small-en-v1.5 (BAAI) — 384-dim |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Web UI | Streamlit (port 8128, NVIDIA black/green theme) |
| REST API | FastAPI + Uvicorn (port 8538) |
| Configuration | Pydantic BaseSettings |
| Testing | pytest (769 tests) |
| Hardware target | NVIDIA DGX Spark (GB10 GPU, 128GB unified, $4,699) |

### 8.2 Service Ports

| Port | Service |
|---|---|
| 8538 | FastAPI REST API |
| 8128 | Streamlit Chat UI |
| 19530 | Milvus vector database (shared) |

### 8.3 Dependencies on HCLS AI Factory

| Dependency | Usage |
|---|---|
| Milvus 2.4 instance | Shared vector database — adds 13 owned collections alongside existing `genomic_evidence` (read-only) |
| `ANTHROPIC_API_KEY` | Loaded from `rag-chat-pipeline/.env` if not set in environment |
| BGE-small-en-v1.5 | Same embedding model as main RAG pipeline |

---

## 9. Demo Scenarios

### 9.1 Validated Demo Queries

**1. "Design an adaptive Phase II protocol for EGFR-mutant NSCLC with osimertinib resistance"**
- Searches: protocols (EGFR NSCLC), endpoints (ORR, PFS), adaptive (biomarker-adaptive), results (AURA3, FLAURA)
- Expected: Biomarker-enriched, seamless Phase II/III, with co-primary ORR/PFS endpoints

**2. "Match patients with BRCA1/2 mutations to open PARP inhibitor trials"**
- Searches: eligibility (BRCA criteria), biomarkers (BRCA1/2), protocols (PARP inhibitor), sites (recruiting)
- Expected: Ranked trial list with eligibility match scores and site proximity

**3. "Identify the top 10 enrolling sites for pediatric ALL trials in the US"**
- Searches: sites (pediatric ALL, US), investigators (pediatric hematology), results (enrollment data)
- Expected: Ranked site list with enrollment rates, investigator profiles, and diversity metrics

**4. "Are there safety signals for hepatotoxicity in the latest checkpoint inhibitor trials?"**
- Searches: safety (hepatotoxicity, checkpoint), literature (hepatic AEs), regulatory (safety reports)
- Expected: Signal detection with frequency, severity, and regulatory response summary

**5. "What is the competitive landscape for GLP-1 receptor agonists in obesity?"**
- Searches: protocols (GLP-1, obesity), results (semaglutide, tirzepatide), regulatory (obesity approvals)
- Expected: Competitive matrix with mechanism, phase, enrollment, and differentiation analysis

---

## 10. File Structure (Actual)

```
clinical_trial_intelligence_agent/
├── src/
│   ├── __init__.py
│   ├── agent.py                     # Autonomous reasoning pipeline
│   ├── models.py                    # Pydantic data models + enums
│   ├── collections.py               # 14 Milvus collection schemas
│   ├── rag_engine.py                # Multi-collection RAG engine
│   ├── clinical_workflows.py        # 10 trial workflows
│   ├── decision_support.py          # 5 decision engines + success estimator
│   ├── knowledge.py                 # Domain knowledge base (40 trials, 13 areas, 9 agencies)
│   ├── query_expansion.py           # 10 expansion maps, 140 aliases
│   ├── cross_modal.py               # Cross-agent integration (4 peer agents)
│   ├── metrics.py                   # Prometheus metrics
│   ├── export.py                    # Report generation (PDF, Markdown, JSON)
│   └── ingest/
│       ├── base.py                  # BaseIngestPipeline
│       ├── protocol_parser.py       # ClinicalTrials.gov REST API v2
│       ├── eligibility_parser.py    # Eligibility criteria extraction
│       ├── safety_parser.py         # Adverse event parser
│       ├── regulatory_parser.py     # FDA/EMA document parser
│       └── literature_parser.py     # PubMed E-utilities
├── app/
│   └── trial_ui.py                 # Streamlit (5 tabs, NVIDIA theme)
├── api/
│   └── main.py                     # FastAPI REST server (26 endpoints)
├── config/
│   └── settings.py                 # Pydantic BaseSettings
├── data/
│   └── reference/                  # Seed data JSON files
├── scripts/
│   ├── setup_collections.py
│   ├── seed_knowledge.py
│   └── validate_e2e.py
├── tests/                          # 769 tests
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```

**46 Python files | ~22,607 lines of code | Apache 2.0**

---

## 11. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Architecture** | Complete | Data models, 14 collection schemas, knowledge base, 5 decision engines, RAG engine, agent orchestrator |
| **Phase 2: Data** | Complete | 40 landmark trials, 13 therapeutic areas, 9 regulatory agencies, 9 adaptive design templates, 140 entity aliases |
| **Phase 3: RAG Integration** | Complete | Multi-collection parallel search, knowledge augmentation, Claude Sonnet 4.6 streaming |
| **Phase 4: Workflows** | Complete | 10 clinical workflows with workflow-specific weight boosting |
| **Phase 5: Testing** | Complete | 769 tests, 100% pass, 0.47s runtime |
| **Phase 6: UI + Demo** | Complete | Streamlit UI on port 8128 (5 tabs), NVIDIA theme, demo scenarios validated |

### Remaining Work

| Item | Priority | Effort |
|---|---|---|
| Real-time ClinicalTrials.gov sync (incremental updates) | Medium | 2-3 days |
| FDA CDER approval database integration | Medium | 1-2 days |
| Multi-region regulatory strategy comparison | Low | 1 week |
| Integration with HCLS AI Factory landing page | Low | 1 hour |

---

## 12. Relationship to HCLS AI Factory

The Clinical Trial Intelligence Agent demonstrates the **translational extension** of the HCLS AI Factory architecture. While the core platform discovers drug candidates (Stage 3), this agent accelerates the path from candidate to clinic by optimizing clinical trial design and execution.

- **Same Milvus instance** — 13 new owned collections alongside existing `genomic_evidence` (read-only)
- **Same embedding model** — BGE-small-en-v1.5 (384-dim)
- **Same LLM** — Claude via Anthropic API
- **Same hardware** — NVIDIA DGX Spark ($4,699)
- **Same patterns** — Pydantic models, BaseIngestPipeline, knowledge graph, query expansion

The platform closes the loop: genomic variants (Stage 1) inform biomarker-enriched trial designs, while trial outcomes feed back into the RAG knowledge base to improve future recommendations.

---

## 13. Credits

- **Adam Jones**
- **Apache 2.0 License**

---

!!! warning "Clinical Decision Support Disclaimer"
    The Clinical Trial Intelligence Agent is a clinical decision support research tool for clinical trial optimization. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals and regulatory experts. Apache 2.0 License.
