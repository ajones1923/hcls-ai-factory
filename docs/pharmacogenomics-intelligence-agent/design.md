# Pharmacogenomics Intelligence Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## 1. Executive Summary

The Pharmacogenomics Intelligence Agent extends the HCLS AI Factory platform to translate patient genetic data into actionable drug prescribing recommendations. It bridges the gap between raw genotype data and clinical prescribing decisions by implementing CPIC-standardized clinical pipelines for star allele calling, phenotype translation, drug-gene matching, phenoconversion detection, HLA hypersensitivity screening, and genotype-guided dosing.

The system searches **15 Milvus vector collections** containing pharmacogene references, CPIC/DPWG clinical guidelines, drug-gene interactions, HLA hypersensitivity associations, phenoconversion models, validated dosing algorithms, clinical evidence, and population allele frequency data. It covers **14 pharmacogenes** interacting with **100+ drugs** across **12 therapeutic categories**, implementing **9 genotype-guided dosing algorithms** and screening **15 HLA-drug hypersensitivity** associations.

The platform enables queries like *"My patient is CYP2D6 *1/*4 and taking codeine — what are the clinical implications?"* that simultaneously search pharmacogene references, CPIC guidelines, dosing algorithms, drug interactions, and population data — returning grounded prescribing recommendations with evidence citations.

### Key Results

| Metric | Value |
|---|---|
| Total Python LOC | **24,577** (17,913 source + 6,664 test) |
| Milvus collections | **15** (14 PGx-specific + 1 shared genomic_evidence) |
| Pharmacogenes in knowledge graph | **14** core (25 extended) |
| Drugs covered | **100+** across 12 therapeutic categories |
| Dosing algorithms | **9** genotype-guided (IWPC warfarin, tacrolimus, fluoropyrimidine, thiopurine, clopidogrel, simvastatin, SSRI, phenytoin, TCA) |
| HLA-drug associations | **15** screened (12 in knowledge graph) |
| CYP inhibitors/inducers | **60** for phenoconversion across CYP2D6, CYP2C19, CYP2C9, CYP3A4/5 |
| Seed data records | **240** across 14 JSON files |
| Clinical workflows | **8** |
| Test suite | **1,001** tests (100% pass, 0.48s runtime) |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | PGx Agent Role |
|---|---|
| **DataStore** | Raw files: CPIC guideline JSONs, PharmVar allele definitions, PharmGKB annotations, FDA PGx labels, population frequency data |
| **DataEngine** | 8 ingest parsers (CPIC, PharmVar, PharmGKB, FDA, population, PubMed, clinical trials, base) |
| **DataBase** | 15 Milvus collections (14 owned + 1 read-only) + 9 knowledge graph dictionaries |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + 6 clinical pipelines + 9 dosing algorithms |
| **AgentEngine** | PGxIntelligenceAgent (Plan-Search-Evaluate-Synthesize) + Streamlit UI (10 tabs) + FastAPI REST |

### 2.2 System Diagram

```
                ┌───────────────────────────────────────────────┐
                │         Streamlit UI (8507)                    │
                │  10 Tabs: Dashboard | Drug Check | Med Review  │
                │  Warfarin | Chemo Safety | HLA Screening |     │
                │  Reports | Evidence | Phenoconversion | PopGen │
                └──────────────────────┬────────────────────────┘
                                       │
                ┌──────────────────────▼────────────────────────┐
                │         FastAPI REST API (8107)                │
                │         16+ endpoints, CORS, Auth              │
                └──────────────────────┬────────────────────────┘
                                       │
                ┌──────────────────────▼────────────────────────┐
                │         PGxIntelligenceAgent                   │
                │  Plan → Search → Evaluate → Synthesize         │
                └──────────────────────┬────────────────────────┘
                                       │
        ┌──────────────────────────────┼──────────────────────────────┐
        │                              │                              │
┌───────▼────────┐          ┌──────────▼──────────┐       ┌──────────▼──────────┐
│ Clinical        │          │ RAG Engine           │       │ Workflow Engine (8) │
│ Pipelines       │          │                      │       │                      │
│                 │          │ BGE-small-en-v1.5    │       │ Drug Check           │
│ Star Allele     │          │ (384-dim embedding)  │       │ Medication Review    │
│ Phenotype Trans.│          │         │            │       │ Warfarin Dosing      │
│ Drug-Gene Match │          │         ▼            │       │ Chemo Safety         │
│ Phenoconversion │          │ Parallel Search      │       │ HLA Screening        │
│ HLA Screening   │          │ 15 Milvus Collections│       │ Phenoconversion      │
│ Genotype Dosing │          │ (ThreadPoolExecutor) │       │ Population Analysis  │
│ (9 algorithms)  │          │         │            │       │ Comprehensive Review │
│                 │          │         ▼            │       │                      │
│                 │          │ Claude Sonnet 4.6    │       │                      │
└───────┬────────┘          └──────────┬──────────┘       └──────────────────────┘
        │                              │
┌───────▼──────────────────────────────▼──────────────────────────────┐
│                  Milvus 2.4 — 15 Collections                        │
│                                                                      │
│  pgx_gene_reference          pgx_drug_guidelines                     │
│  pgx_drug_interactions       pgx_hla_hypersensitivity                │
│  pgx_phenoconversion         pgx_dosing_algorithms                   │
│  pgx_clinical_evidence       pgx_population_frequencies              │
│  pgx_clinical_trials         pgx_fda_labels                          │
│  pgx_therapeutic_alternatives pgx_patient_diplotypes                 │
│  pgx_implementation          pgx_education                           │
│  genomic_evidence [shared, read-only]                                │
└──────────────────────────────────────────────────────────────────────┘
```

---

## 3. Data Collections — Actual State

### 3.1 Collection Catalog

| # | Collection | Records | Weight | Primary Use |
|---|---|---|---|---|
| 1 | `pgx_gene_reference` | 25 | 0.10 | Pharmacogene definitions, allele tables |
| 2 | `pgx_drug_guidelines` | 35 | 0.14 | CPIC/DPWG clinical guidelines |
| 3 | `pgx_drug_interactions` | 30 | 0.12 | Drug-gene interaction severity |
| 4 | `pgx_hla_hypersensitivity` | 15 | 0.10 | HLA allele-drug reaction screening |
| 5 | `pgx_phenoconversion` | 20 | 0.08 | CYP inhibitor/inducer phenotype shifts |
| 6 | `pgx_dosing_algorithms` | 9 | 0.07 | Genotype-guided dose calculations |
| 7 | `pgx_clinical_evidence` | 25 | 0.06 | Published pharmacogenomic studies |
| 8 | `pgx_population_frequencies` | 18 | 0.05 | Population allele frequency data |
| 9 | `pgx_clinical_trials` | 15 | 0.05 | PGx-relevant clinical trials |
| 10 | `pgx_fda_labels` | 15 | 0.06 | FDA pharmacogenomic label annotations |
| 11 | `pgx_therapeutic_alternatives` | 12 | 0.05 | Alternative drug recommendations |
| 12 | `pgx_patient_diplotypes` | 8 | 0.04 | Sample patient diplotype profiles |
| 13 | `pgx_implementation` | 8 | 0.03 | Clinical implementation protocols |
| 14 | `pgx_education` | 5 | 0.02 | Educational resources |
| 15 | `genomic_evidence` | shared | 0.03 | Cross-modal genomic queries |
| | **Total seed records** | **240** | **1.00** | |

### 3.2 Index Configuration

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 128 (all collections — small record counts) |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Clinical Pipelines

### 4.1 Star Allele Calling

Translates VCF-derived pharmacogene variants into star allele nomenclature following PharmVar conventions. Supports diplotype assignment for all 14 core pharmacogenes.

### 4.2 Phenotype Translation

Maps diplotypes to CPIC-standardized metabolizer phenotypes using activity score methodology:

| Gene | Activity Score System | Example |
|---|---|---|
| CYP2D6 | AS 0-2+ scale | *1/*4 → AS 1.0 → Intermediate Metabolizer |
| CYP2C19 | AS 0-2 scale | *1/*2 → AS 1.0 → Intermediate Metabolizer |
| CYP2C9 | AS 0-2 scale | *1/*3 → AS 1.0 → Intermediate Metabolizer |
| DPYD | AS 0-2 scale | *1/*2A → AS 1.0 → Intermediate Metabolizer |
| TPMT | AS 0-2 scale | *1/*3A → AS 1.0 → Intermediate Metabolizer |

### 4.3 Drug-Gene Matching

CPIC guideline lookup with alert severity classification:

| Alert Level | Meaning | Example |
|---|---|---|
| **Contraindicated** | Do not prescribe | CYP2D6 Ultra-Rapid + codeine |
| **Major** | Significant dose adjustment needed | CYP2C19 Poor Metabolizer + clopidogrel |
| **Moderate** | Consider alternative or monitor | CYP2C9 Intermediate + warfarin |
| **Informative** | No action needed | Normal metabolizer + standard drug |

### 4.4 Phenoconversion Detection

Models CYP inhibitor/inducer effects on metabolizer phenotype:

| CYP Enzyme | Strong Inhibitors | Strong Inducers |
|---|---|---|
| CYP2D6 | Fluoxetine, paroxetine, quinidine, bupropion | None clinically significant |
| CYP2C19 | Omeprazole, esomeprazole, fluvoxamine | Rifampin, St. John's wort |
| CYP2C9 | Fluconazole, amiodarone | Rifampin, phenytoin |
| CYP3A4/5 | Ketoconazole, itraconazole, clarithromycin | Rifampin, carbamazepine, phenytoin |

### 4.5 HLA Screening

Pre-prescription screening for 15 HLA-drug hypersensitivity associations:

| HLA Allele | Drug | Reaction | Prevalence |
|---|---|---|---|
| HLA-B*57:01 | Abacavir | Hypersensitivity syndrome | 5-8% Caucasian |
| HLA-B*58:01 | Allopurinol | SJS/TEN | 6-8% Southeast Asian |
| HLA-B*15:02 | Carbamazepine | SJS/TEN | 8-15% Southeast Asian |
| HLA-A*31:01 | Carbamazepine | DRESS/maculopapular | 2-5% Caucasian |
| HLA-B*57:01 | Flucloxacillin | DILI | 5-8% Caucasian |

### 4.6 Genotype-Guided Dosing Algorithms (9)

| Algorithm | Input | Output | Validation Source |
|---|---|---|---|
| IWPC Warfarin | CYP2C9, VKORC1, age, weight, race | Weekly dose (mg) | IWPC 2009 |
| CYP3A5 Tacrolimus | CYP3A5 genotype, weight | Starting dose (mg/day) | CPIC 2015 |
| DPYD Fluoropyrimidine | DPYD genotype, BSA | Dose reduction (%) | CPIC 2017 |
| TPMT/NUDT15 Thiopurine | TPMT/NUDT15 genotype | Starting dose (mg/m2) | CPIC 2018 |
| CYP2C19 Clopidogrel | CYP2C19 genotype | Alternative selection | CPIC 2013 |
| SLCO1B1 Simvastatin | SLCO1B1 genotype | Max dose / alternative | CPIC 2014 |
| CYP2D6/2C19 SSRI | CYP2D6/CYP2C19 genotype | Dose adjustment (%) | CPIC 2015 |
| CYP2C9 Phenytoin | CYP2C9, HLA-B*15:02 | Dose + safety screen | CPIC 2014 |
| CYP2D6 TCA | CYP2D6 genotype | Dose adjustment (%) | CPIC 2016 |

---

## 5. Multi-Collection RAG Engine

### 5.1 Search Flow

```
User Query: "CYP2D6 *1/*4 patient on codeine — recommendations?"
    │
    ├── 1. Embed query with BGE asymmetric prefix               [< 5 ms]
    │
    ├── 2. Parallel search across 15 collections (top-5 each)   [10-16 ms]
    │   ├── pgx_drug_guidelines:   CYP2D6 codeine CPIC          (score: 0.85-0.92)
    │   ├── pgx_gene_reference:    CYP2D6 allele definitions    (score: 0.80-0.88)
    │   ├── pgx_drug_interactions: CYP2D6 opioid interactions   (score: 0.78-0.86)
    │   ├── pgx_dosing_algorithms: CYP2D6 dose adjustment       (score: 0.72-0.80)
    │   └── pgx_therapeutic_alternatives: Non-CYP2D6 analgesics (score: 0.70-0.78)
    │
    ├── 3. Clinical pipeline: Star allele → phenotype            [< 10 ms]
    │      *1/*4 → AS 1.0 → Intermediate Metabolizer
    │
    ├── 4. Drug-gene matching                                    [< 5 ms]
    │      CYP2D6 IM + codeine → Reduced efficacy (less morphine conversion)
    │      Alert: MODERATE — consider tramadol alternative with dose adjustment
    │
    ├── 5. Knowledge base augmentation                           [< 1 ms]
    │
    └── 6. Stream Claude Sonnet 4.6 response                    [~20-24 sec]
           CPIC-grounded recommendation with alternative drugs,
           dose adjustments, and monitoring plan
```

**Total: ~24 sec** (dominated by LLM generation; retrieval + clinical pipelines ~30 ms)

---

## 6. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 6.1 Clinical Pipeline Performance

| Pipeline | Latency |
|---|---|
| Star allele calling (per gene) | <5 ms |
| Phenotype translation (all 14 genes) | <10 ms |
| Drug-gene matching (full medication list) | <20 ms |
| Phenoconversion detection (60 inhibitors) | <15 ms |
| HLA screening (15 associations) | <10 ms |
| Dosing algorithm (single drug) | <10 ms |
| **All pipelines combined** | **<70 ms** |

### 6.2 RAG Query Performance

| Operation | Latency |
|---|---|
| Full query (retrieve + Claude generate) | ~24 sec |
| Streaming query (time to first token) | ~3 sec |
| 15-collection parallel search | 10-16 ms |
| Response length | 800-2000 chars |

### 6.3 Seed Performance

| Operation | Duration |
|---|---|
| Seed all 15 collections (240 vectors) | ~1.5 min |
| BGE-small embedding per batch | ~1.5 sec |
| Milvus insert per collection | <500 ms |

---

## 7. Infrastructure

### 7.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4, localhost:19530 |
| Embeddings | BGE-small-en-v1.5 (BAAI) — 384-dim |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Web UI | Streamlit (port 8507, NVIDIA black/green theme) |
| REST API | FastAPI + Uvicorn (port 8107) |
| Configuration | Pydantic BaseSettings |
| Testing | pytest (1,001 tests) |
| Hardware target | NVIDIA DGX Spark (GB10 GPU, 128GB unified, $4,699) |

### 7.2 Service Ports

| Port | Service |
|---|---|
| 8107 | FastAPI REST API |
| 8507 | Streamlit Chat UI |
| 19530 | Milvus vector database (shared) |

### 7.3 Dependencies on HCLS AI Factory

| Dependency | Usage |
|---|---|
| Milvus 2.4 instance | Shared vector database — adds 14 owned collections alongside existing `genomic_evidence` (read-only) |
| `ANTHROPIC_API_KEY` | Shared Anthropic API key |
| BGE-small-en-v1.5 | Same embedding model as main RAG pipeline |

---

## 8. Demo Scenarios

### 8.1 Validated Demo Queries

**1. "CYP2D6 *1/*4 patient prescribed codeine for post-surgical pain"**
- Phenotype: Intermediate Metabolizer (AS 1.0)
- Alert: Reduced codeine → morphine conversion, consider alternative analgesic
- Dosing: Standard dose may be ineffective; tramadol or non-opioid alternative recommended

**2. "Screen HLA-B*57:01 before starting abacavir for HIV"**
- HLA Screener: HLA-B*57:01 positive → Abacavir CONTRAINDICATED
- Alternative: Tenofovir-based regimen recommended

**3. "Calculate warfarin dose for CYP2C9 *1/*3, VKORC1 AG, 70kg, 65yo, Caucasian male"**
- IWPC Algorithm: Estimated weekly dose with CYP2C9 intermediate sensitivity adjustment
- VKORC1 AG: Intermediate warfarin sensitivity

**4. "Review medication list for PGx interactions: omeprazole, clopidogrel, metoprolol, sertraline"**
- Phenoconversion: Omeprazole inhibits CYP2C19 → affects clopidogrel activation
- Drug-gene: CYP2D6 status affects metoprolol and sertraline metabolism
- Alert: Omeprazole + clopidogrel interaction flagged

**5. "DPYD *1/*2A patient starting 5-FU chemotherapy — dose adjustment?"**
- Phenotype: Intermediate Metabolizer (AS 1.0)
- CPIC: 50% dose reduction for fluoropyrimidines (CONTRAINDICATED if Poor Metabolizer)
- Monitoring: DPD enzyme activity test recommended for confirmation

---

## 9. File Structure (Actual)

```
pharmacogenomics_intelligence_agent/
├── src/
│   ├── agent.py                     # PGxIntelligenceAgent orchestrator
│   ├── models.py                    # Pydantic data models + enums
│   ├── collections.py               # 15 Milvus collection schemas
│   ├── rag_engine.py                # PGxRAGEngine multi-collection search
│   ├── clinical_workflows.py        # 8 clinical workflows
│   ├── clinical_pipelines.py        # Star allele, phenotype, dosing pipelines
│   ├── hla_screener.py              # 15 HLA-drug screening associations
│   ├── phenoconversion.py           # CYP inhibitor/inducer phenotype modeling
│   ├── dosing_algorithms.py         # 9 genotype-guided dosing calculators
│   ├── knowledge.py                 # 9 knowledge graph dictionaries
│   ├── query_expansion.py           # PGx synonym expansion
│   ├── cross_modal.py               # Cross-agent integration
│   ├── metrics.py                   # Prometheus metrics
│   ├── export.py                    # Report generation
│   └── ingest/
│       ├── cpic_parser.py           # CPIC guideline ingest
│       ├── pharmvar_parser.py       # PharmVar allele definitions
│       ├── pharmgkb_parser.py       # PharmGKB annotations
│       ├── fda_parser.py            # FDA PGx label annotations
│       ├── population_parser.py     # Population frequency data
│       ├── pubmed_parser.py         # PubMed literature
│       └── trials_parser.py         # Clinical trial data
├── app/
│   └── pgx_ui.py                   # Streamlit (10 tabs, NVIDIA theme)
├── api/
│   └── main.py                     # FastAPI REST server (16+ endpoints)
├── config/
│   └── settings.py                 # Pydantic BaseSettings
├── data/
│   └── reference/                  # 14 JSON seed files (240 records)
├── tests/                          # 1,001 tests across 16 test files
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```

**83 files | ~24,577 lines of code | Apache 2.0**

---

## 10. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Architecture** | Complete | 15 collections, knowledge graph, 6 clinical pipelines, 9 dosing algorithms, 8 workflows |
| **Phase 2: Data** | Complete | 240 seed records across 14 JSON files, 14 pharmacogenes, 100+ drugs, 15 HLA associations |
| **Phase 3: RAG Integration** | Complete | Multi-collection parallel search, Claude Sonnet 4.6 streaming |
| **Phase 4: Testing** | Complete | 1,001 tests, 100% pass, 0.48s runtime |
| **Phase 5: UI + Demo** | Complete | 10-tab Streamlit UI, NVIDIA theme, 5 demo scenarios validated |

### Remaining Work

| Item | Priority | Effort |
|---|---|---|
| CPIC guideline auto-update pipeline | Medium | 2-3 days |
| Additional population frequency data (gnomAD v4) | Medium | 1-2 days |
| EHR integration via FHIR PGx profiles | Low | 1 week |
| Integration with HCLS AI Factory landing page | Low | 1 hour |

---

## 11. Relationship to HCLS AI Factory

The Pharmacogenomics Intelligence Agent demonstrates the **prescribing translation** capability of the HCLS AI Factory architecture. While Stage 3 discovers drug candidates, this agent ensures those drugs are prescribed safely and effectively based on individual patient genetics.

- **Same Milvus instance** — 14 new owned collections alongside existing `genomic_evidence` (read-only)
- **Same embedding model** — BGE-small-en-v1.5 (384-dim)
- **Same LLM** — Claude via Anthropic API
- **Same hardware** — NVIDIA DGX Spark ($4,699)
- **Same patterns** — Pydantic models, BaseIngestPipeline, knowledge graph, query expansion

The critical integration point: Stage 1 (genomics) produces VCF files containing pharmacogene variants; this agent translates those variants into prescribing guidance, closing the loop from DNA to dosing decision.

---

## 12. Credits

- **Adam Jones**
- **Apache 2.0 License**

---

!!! warning "Clinical Decision Support Disclaimer"
    The Pharmacogenomics Intelligence Agent is a clinical decision support research tool for pharmacogenomic interpretation. It is not FDA-cleared and is not intended as a standalone diagnostic device. All prescribing recommendations should be reviewed by qualified healthcare professionals and pharmacists. CPIC guideline compliance should be verified against current published guidelines. Apache 2.0 License.
