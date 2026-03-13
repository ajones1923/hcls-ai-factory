# Pharmacogenomics Intelligence Agent -- Project Bible

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0
**Repository:** hcls-ai-factory/ai_agent_adds/pharmacogenomics_intelligence_agent

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Vision & Mission](#2-vision--mission)
3. [System Overview](#3-system-overview)
4. [File Inventory](#4-file-inventory)
5. [Collections Catalog](#5-collections-catalog)
6. [Knowledge Graph](#6-knowledge-graph)
7. [Clinical Pipelines](#7-clinical-pipelines)
8. [Clinical Workflows](#8-clinical-workflows)
9. [RAG Engine](#9-rag-engine)
10. [Agent Architecture](#10-agent-architecture)
11. [Query Expansion](#11-query-expansion)
12. [Data Models](#12-data-models)
13. [Export System](#13-export-system)
14. [API Reference](#14-api-reference)
15. [UI Guide](#15-ui-guide)
16. [Metrics & Monitoring](#16-metrics--monitoring)
17. [Scheduler](#17-scheduler)
18. [Configuration](#18-configuration)
19. [Docker Deployment](#19-docker-deployment)
20. [Testing](#20-testing)
21. [Port Map](#21-port-map)
22. [Tech Stack](#22-tech-stack)
23. [Future Roadmap](#23-future-roadmap)

---

## 1. Executive Summary

The **Pharmacogenomics Intelligence Agent** is a domain-specialized retrieval-augmented generation (RAG) system that translates patient genetic data into actionable drug prescribing recommendations. It searches 15 Milvus vector collections containing pharmacogene references, CPIC/DPWG clinical guidelines, drug-gene interactions, HLA hypersensitivity associations, phenoconversion models, validated dosing algorithms, clinical evidence, population allele frequency data, clinical trial records, FDA pharmacogenomic labels, therapeutic alternatives, patient diplotype profiles, implementation protocols, and educational resources.

The system implements nine genotype-guided dosing algorithms (IWPC warfarin, CYP3A5 tacrolimus, DPYD fluoropyrimidine, TPMT+NUDT15 thiopurine, CYP2C19 clopidogrel, SLCO1B1 simvastatin, CYP2D6/CYP2C19 SSRI, CYP2C9 phenytoin, CYP2D6 TCA), screens 15 HLA-drug hypersensitivity associations, models phenoconversion across 30+ CYP inhibitors/inducers, and covers 25 pharmacogenes interacting with 100+ drugs across 12 therapeutic categories.

The Pharmacogenomics Intelligence Agent is built as part of the HCLS AI Factory, a three-stage precision medicine platform (Genomics, RAG/Chat, Drug Discovery) designed to run end-to-end on a single NVIDIA DGX Spark.

### Codebase at a Glance

| Metric | Value |
|--------|-------|
| Total Python LOC | 24,577 (17,913 source + 6,664 test) |
| Total files | 83 (53 Python, 14 JSON seed, 8 config/infra, 8 docs) |
| Milvus collections | 15 (14 PGx-specific + 1 shared genomic_evidence) |
| Seed data records | 240 across 14 JSON files |
| Knowledge graph pharmacogenes | 25 |
| Drugs covered | 100+ across 12 therapeutic categories |
| HLA-drug associations | 15 (screener) / 12 (knowledge graph) |
| CYP inhibitors/inducers | 30+ for phenoconversion |
| Dosing algorithms | 9 genotype-guided |
| Clinical workflows | 8 |
| Test suite | 1,001 tests, all passing in 0.48s |
| LLM | Claude Sonnet 4.6 (claude-sonnet-4-6) |
| Embedding model | BGE-small-en-v1.5 (384-dim) |

---

## 2. Vision & Mission

### The Problem

Adverse drug reactions (ADRs) affect approximately 95% of the population through genetic variation in drug-metabolizing enzymes, transporters, and HLA alleles. ADR-related costs exceed $528 billion per year globally, including hospitalizations, treatment failures, and preventable deaths. Despite the availability of CPIC guidelines for over 100 gene-drug pairs, fewer than 5% of patients receive pharmacogenomic testing before being prescribed high-risk medications.

### Mission Statement

Democratize pharmacogenomic decision support by providing unified, evidence-grounded intelligence across the entire PGx landscape -- from star allele interpretation through genotype-guided dosing -- on hardware that any clinic, pharmacy, or research lab can afford.

### Design Principles

1. **Clinically Validated**: Every recommendation traces to CPIC, DPWG, or FDA guidelines with evidence level classification.
2. **Safety-First**: Contraindicated combinations generate immediate alerts; HLA screening flags are never suppressed.
3. **Population-Aware**: Allele frequency data spans all major populations; equity implications are surfaced.
4. **Phenoconversion-Aware**: Concomitant medications are checked for CYP inhibition/induction before phenotype assignment.
5. **Quantitative**: Dosing algorithms produce specific dose recommendations, not vague guidance.

---

## 3. System Overview

### High-Level Architecture

```
+------------------------------------------------------------------+
|                    STREAMLIT UI (8507)                             |
|  [Dashboard] [Drug Check] [Med Review] [Warfarin] [Chemo Safety] |
|  [HLA Screening] [Reports] [Evidence] [Phenoconversion] [PopGen] |
|  Sidebar: Gene filter, Drug filter, Collection stats, Export      |
+------------------------------------------------------------------+
         |                                          |
         v                                          v
+------------------+                    +---------------------+
| PGxIntelligence  |                    | FastAPI REST (8107) |
|     Agent        |                    | 16+ endpoints       |
|  plan -> search  |                    | /health /collections|
|  -> evaluate     |                    | /query /search      |
|  -> synthesize   |                    | /v1/pgx/* (7 PGx)  |
+------------------+                    +---------------------+
         |                                          |
         v                                          v
+------------------------------------------------------------------+
|                      PGxRAGEngine                                 |
|  embed -> parallel search -> merge & rank -> knowledge augment    |
|  -> clinical pipeline (star allele, phenoconversion, HLA, dosing) |
|  -> LLM synthesis (Claude Sonnet 4.6)                             |
+------------------------------------------------------------------+
         |
         v
+------------------------------------------------------------------+
|                 PGxCollectionManager                              |
|  15 collections | IVF_FLAT index | COSINE metric                 |
|  ThreadPoolExecutor parallel search across all collections        |
+------------------------------------------------------------------+
         |
         v
+------------------------------------------------------------------+
|                    Milvus 2.4 (19530)                             |
|  BGE-small-en-v1.5 embeddings (384-dim)                          |
|                                                                   |
|  pgx_gene_reference      pgx_drug_guidelines                     |
|  pgx_drug_interactions   pgx_hla_hypersensitivity                |
|  pgx_phenoconversion     pgx_dosing_algorithms                   |
|  pgx_clinical_evidence   pgx_population_data                     |
|  pgx_clinical_trials     pgx_fda_labels                          |
|  pgx_drug_alternatives   pgx_patient_profiles                    |
|  pgx_implementation      pgx_education                           |
|  genomic_evidence (read-only, shared)                             |
+------------------------------------------------------------------+
```

---

## 4. File Inventory

### Source Files (src/)

| File | Lines | Purpose |
|------|-------|---------|
| `knowledge.py` | 2,657 | 25 pharmacogenes, 12 drug categories, 12 HLA associations, CYP inhibitors/inducers, drug alternatives, activity score tables, entity aliases |
| `collections.py` | 1,547 | 15 Milvus collection schemas, parallel search, IVF_FLAT indexing |
| `export.py` | 1,307 | Markdown, JSON, PDF, FHIR R4 export with PGx alerts and drug interaction matrix |
| `query_expansion.py` | 1,254 | 14 domain-specific expansion maps for drug, gene, phenotype, HLA, dosing terms |
| `pgx_pipeline.py` | 1,600 | StarAlleleCaller, PhenotypeTranslator, DrugGeneMatcher -- VCF-to-PGx pipeline |
| `rag_engine.py` | 799 | Multi-collection RAG with PGx knowledge augmentation and clinical alerts |
| `agent.py` | 588 | Autonomous plan-search-evaluate-synthesize-report reasoning agent |
| `dosing.py` | 1,499 | 9 genotype-guided dosing algorithms (IWPC warfarin, tacrolimus, fluoropyrimidine, thiopurine, clopidogrel, simvastatin, SSRI, phenytoin, TCA) |
| `hla_screener.py` | 725 | 15 HLA-drug hypersensitivity screening with population-specific risk |
| `models.py` | 616 | Enums, collection models, query/response schemas, alert types |
| `phenoconversion.py` | 517 | CYP inhibitor/inducer phenoconversion detection and adjustment |
| `metrics.py` | 399 | 22 Prometheus metrics (10 histograms, 8 counters, 4 gauges) |
| `scheduler.py` | 232 | Weekly PubMed + ClinicalTrials.gov automated ingest refresh |

### Ingest Parsers (src/ingest/)

| File | Purpose |
|------|---------|
| `base.py` | Abstract base class for all ingest parsers |
| `cpic_parser.py` | CPIC guideline API ingest |
| `pharmvar_parser.py` | PharmVar star allele database ingest |
| `pharmgkb_parser.py` | PharmGKB clinical annotation ingest |
| `fda_label_parser.py` | FDA pharmacogenomic drug label ingest |
| `population_parser.py` | Population allele frequency data ingest |
| `pubmed_parser.py` | PubMed PGx literature ingest |
| `clinical_trials_parser.py` | ClinicalTrials.gov PGx trial ingest |

### API Layer (api/)

| File | Lines | Purpose |
|------|-------|---------|
| `main.py` | 628 | FastAPI app with 8 core endpoints, CORS, lifespan management |
| `routes/pgx_clinical.py` | 858 | 7 PGx-specific clinical decision support endpoints |
| `routes/reports.py` | — | Report generation endpoints |
| `routes/events.py` | — | Event audit trail endpoints |

### UI Layer (app/)

| File | Lines | Purpose |
|------|-------|---------|
| `pgx_ui.py` | 2,152 | 10-tab Streamlit interface with NVIDIA dark theme |

### Configuration (config/)

| File | Lines | Purpose |
|------|-------|---------|
| `settings.py` | 191 | Pydantic BaseSettings with PGX_ env_prefix, 15 collection weights |

### Seed Data (data/reference/)

14 JSON files with 240 total records:

| File | Content |
|------|---------|
| `gene_reference.json` | Pharmacogene star allele definitions |
| `drug_guidelines.json` | CPIC/DPWG prescribing guidelines |
| `drug_interactions.json` | Drug-gene interaction records |
| `hla_hypersensitivity.json` | HLA-drug ADR associations |
| `phenoconversion.json` | CYP inhibitor/inducer phenoconversion data |
| `dosing_algorithms.json` | Genotype-guided dosing formulas |
| `clinical_evidence.json` | Published PGx outcome studies |
| `population_data.json` | Population allele frequencies |
| `clinical_trials.json` | PGx clinical trial records |
| `fda_labels.json` | FDA pharmacogenomic labeling |
| `drug_alternatives.json` | Genotype-guided therapeutic alternatives |
| `patient_profiles.json` | Patient diplotype-phenotype profiles |
| `implementation.json` | Clinical PGx implementation programs |
| `education.json` | PGx educational resources |

### Tests (tests/)

| File | Purpose |
|------|---------|
| `test_knowledge.py` | Knowledge graph entries and entity resolution |
| `test_collections.py` | Collection schemas and parallel search |
| `test_models.py` | Pydantic model validation |
| `test_rag_engine.py` | RAG retrieval and synthesis |
| `test_agent.py` | Agent reasoning pipeline |
| `test_pgx_pipeline.py` | Star allele calling, phenotype translation, drug matching |
| `test_dosing.py` | Dosing algorithm calculations |
| `test_hla_screener.py` | HLA screening logic |
| `test_phenoconversion.py` | Phenoconversion detection |
| `test_query_expansion.py` | Query expansion maps |
| `test_export.py` | Export rendering (Markdown, JSON, PDF, FHIR) |
| `test_metrics.py` | Prometheus metric instrumentation |
| `test_settings.py` | Configuration loading |
| `test_scheduler.py` | Scheduler job registration |
| `test_api.py` | API endpoint testing |
| `test_ingest.py` | Ingest parser testing |
| `conftest.py` | Shared fixtures for all test modules |

### Infrastructure

| File | Purpose |
|------|---------|
| `Dockerfile` | Python 3.12 slim image |
| `docker-compose.yml` | 6-service stack (etcd, MinIO, Milvus, Streamlit, API, setup) |
| `requirements.txt` | Python dependencies |
| `README.md` | Project overview |

---

## 5. Collections Catalog

### 5.1 Collection Summary

| # | Collection | Weight | Description | Key Fields |
|---|-----------|--------|-------------|------------|
| 1 | `pgx_gene_reference` | 0.10 | Pharmacogene star allele definitions & activity scores | gene, star_allele, defining_variants, activity_score, function_status, allele_frequency_global |
| 2 | `pgx_drug_guidelines` | 0.14 | CPIC/DPWG clinical prescribing guidelines | gene, drug, guideline_body, cpic_level, phenotype, recommendation |
| 3 | `pgx_drug_interactions` | 0.12 | Drug-gene interaction records from PharmGKB | gene, drug, interaction_type, evidence_level, pk_effect |
| 4 | `pgx_hla_hypersensitivity` | 0.10 | HLA-mediated adverse drug reaction screening | hla_allele, drug, reaction_type, severity, population_risk |
| 5 | `pgx_phenoconversion` | 0.08 | Metabolic phenoconversion via drug-drug interactions | inhibitor_drug, affected_enzyme, inhibitor_strength, phenotype_shift |
| 6 | `pgx_dosing_algorithms` | 0.07 | Genotype-guided dosing algorithms & formulas | drug, algorithm_name, gene, dose_adjustment, formula |
| 7 | `pgx_clinical_evidence` | 0.08 | Published PGx clinical evidence & outcomes | pmid, study_type, gene, drug, outcome, evidence_level |
| 8 | `pgx_population_data` | 0.06 | Population-specific allele frequency data | gene, allele, population, frequency, sample_size |
| 9 | `pgx_clinical_trials` | 0.04 | PGx-related clinical trials | nct_id, gene, drug, phase, status, enrollment |
| 10 | `pgx_fda_labels` | 0.06 | FDA pharmacogenomic labeling information | drug, gene, label_section, biomarker_status, fda_action |
| 11 | `pgx_drug_alternatives` | 0.05 | Genotype-guided therapeutic alternatives | gene, phenotype, drug_to_avoid, alternative_drug, rationale |
| 12 | `pgx_patient_profiles` | 0.03 | Patient diplotype-phenotype profiles | patient_id, gene, diplotype, phenotype, activity_score |
| 13 | `pgx_implementation` | 0.02 | Clinical PGx implementation programs | institution, program_type, genes_tested, ehr_integration |
| 14 | `pgx_education` | 0.02 | PGx educational resources & guidelines | topic, target_audience, content_type, learning_objectives |
| 15 | `genomic_evidence` | 0.03 | Shared read-only genomic variants (from rag-chat-pipeline) | chrom, pos, ref, alt, gene, consequence |

### 5.2 Index Configuration

- **Index type**: IVF_FLAT
- **Distance metric**: COSINE
- **nlist**: 1024
- **Embedding dimension**: 384 (BGE-small-en-v1.5)

---

## 6. Knowledge Graph

**File:** `src/knowledge.py` (2,657 lines)

### 6.1 Nine Structured Dictionaries

| Dictionary | Entries | Description |
|-----------|---------|-------------|
| `PHARMACOGENES` | 25 | Gene entries with clinical pharmacogenomic data (substrates, variants, CPIC guidelines) |
| `METABOLIZER_PHENOTYPES` | 5 | Phenotype classifications (UM, RM, NM, IM, PM) with activity score ranges |
| `DRUG_CATEGORIES` | 12 | Therapeutic categories with member drugs (opioids, antidepressants, anticoagulants, etc.) |
| `CYP_INHIBITORS` | 4 enzymes | Strong/moderate/weak inhibitors per CYP enzyme |
| `CYP_INDUCERS` | 3 enzymes | Strong/moderate inducers per CYP enzyme |
| `HLA_DRUG_ASSOCIATIONS` | 12 | HLA-drug hypersensitivity pairs with reaction type and severity |
| `DRUG_ALTERNATIVES` | 60 | Gene-phenotype-specific drug substitutions with rationale |
| `ACTIVITY_SCORE_TABLES` | 8 genes | CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, UGT1A1 activity score-to-phenotype mappings |
| `ENTITY_ALIASES` | 116 | Aliases for drugs (brand/generic), genes, and phenotypes |

### 6.2 The 25 Pharmacogenes

**CYP Enzymes (8):** CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP3A5, CYP2B6, CYP1A2, CYP4F2

**Phase II Enzymes (4):** UGT1A1, NAT2, TPMT, DPYD

**Transporters (2):** SLCO1B1, ABCB1

**HLA Genes (3):** HLA-A, HLA-B, HLA-DRB1

**Other Pharmacogenes (8):** VKORC1, NUDT15, G6PD, IFNL3, RYR1, CACNA1S, CYP2A6, COMT

### 6.3 The 12 Therapeutic Drug Categories

1. Opioids (codeine, tramadol, oxycodone, hydrocodone)
2. Antidepressants (escitalopram, citalopram, sertraline, amitriptyline, nortriptyline, venlafaxine)
3. Anticoagulants (warfarin)
4. Antiplatelets (clopidogrel, prasugrel, ticagrelor)
5. Antiepileptics (carbamazepine, phenytoin, oxcarbazepine)
6. Antivirals (abacavir, efavirenz)
7. Immunosuppressants (tacrolimus, azathioprine, mercaptopurine)
8. Statins (simvastatin, atorvastatin, rosuvastatin)
9. Proton Pump Inhibitors (omeprazole, lansoprazole, pantoprazole)
10. Fluoropyrimidines (fluorouracil, capecitabine)
11. Thiopurines (azathioprine, mercaptopurine, thioguanine)
12. Antipsychotics (aripiprazole, risperidone, haloperidol)

---

## 7. Clinical Pipelines

### 7.1 Star Allele Calling (`pgx_pipeline.py` -- StarAlleleCaller)

Resolves VCF variant calls to pharmacogene star allele nomenclature using PharmVar-aligned definitions. Handles:
- Simple SNV-based allele calls (e.g., CYP2C19*2 from rs4244285)
- Compound heterozygous diplotype resolution
- Gene deletion detection (CYP2D6*5)
- Gene duplication detection (CYP2D6*1xN, *2xN)

### 7.2 Phenotype Translation (`pgx_pipeline.py` -- PhenotypeTranslator)

Converts star allele diplotypes to CPIC-standardized metabolizer phenotypes using activity score summation:
- Ultra-Rapid Metabolizer (UM): Activity score > 2.0
- Rapid Metabolizer (RM): Activity score 1.5-2.0
- Normal Metabolizer (NM): Activity score 1.0-2.0
- Intermediate Metabolizer (IM): Activity score 0.5-1.0
- Poor Metabolizer (PM): Activity score 0.0

### 7.3 Drug-Gene Matching (`pgx_pipeline.py` -- DrugGeneMatcher)

Cross-references patient metabolizer profiles against active medication lists:
- CPIC Level A/B guideline lookup per gene-drug pair
- Alert severity classification (Contraindicated, Major, Moderate, Minor, Informational)
- Therapeutic alternative identification with evidence level

### 7.4 Phenoconversion Detection (`phenoconversion.py`)

Adjusts genetic phenotype based on concomitant CYP inhibitors/inducers:
- 30+ inhibitors across CYP2D6, CYP2C19, CYP2C9, CYP3A4/5
- FDA-classified inhibition strength (strong, moderate, weak)
- Phenotype downshift modeling (e.g., NM + strong inhibitor = PM)

### 7.5 HLA Screening (`hla_screener.py`)

Pre-emptive screening for 15 HLA-drug hypersensitivity associations:
- HLA-B*57:01 / abacavir (hypersensitivity syndrome)
- HLA-B*15:02 / carbamazepine (SJS/TEN)
- HLA-B*15:02 / phenytoin (SJS/TEN)
- HLA-A*31:01 / carbamazepine (DRESS/SJS)
- Population-specific risk stratification (e.g., Southeast Asian prevalence)

### 7.6 Dosing Algorithms (`dosing.py`)

Nine genotype-guided dosing calculators:
1. **IWPC Warfarin** -- Age, height, weight, race, CYP2C9, VKORC1, amiodarone, enzyme inducer status
2. **CYP3A5 Tacrolimus** -- Expresser (*1 carrier) vs non-expresser (*3/*3) starting dose
3. **DPYD Fluoropyrimidine** -- Activity score-based dose reduction (0%, 25%, 50%, 100% reduction)
4. **TPMT+NUDT15 Thiopurine** -- Combined activity score for azathioprine/mercaptopurine dosing
5. **CYP2C19 Clopidogrel** -- Metabolizer status-based alternative antiplatelet selection
6. **SLCO1B1 Simvastatin** -- Transporter function-guided statin dose cap and alternative selection
7. **CYP2D6/CYP2C19 SSRI** -- Dual-gene antidepressant dosing with metabolizer-specific adjustments
8. **CYP2C9 Phenytoin** -- Genotype-guided anticonvulsant dose reduction for decreased-function alleles
9. **CYP2D6 TCA** -- Tricyclic antidepressant dosing based on CYP2D6 metabolizer status

---

## 8. Clinical Workflows

The system supports 8 clinical workflows accessible via the UI and API:

| # | Workflow | Tab / Endpoint | Description |
|---|----------|---------------|-------------|
| 1 | Pre-emptive PGx Panel | Dashboard | Multi-gene panel interpretation with activity scores and phenotypes |
| 2 | Opioid Safety Review | Drug Check | CYP2D6-guided opioid prescribing with prodrug metabolism assessment |
| 3 | Anticoagulant Optimization | Warfarin Dosing | IWPC algorithm with CYP2C9/VKORC1 genotype-guided dose calculation |
| 4 | Antidepressant Selection | Drug Check | CYP2D6/CYP2C19-guided SSRI/SNRI/TCA selection |
| 5 | Statin Myopathy Prevention | Drug Check | SLCO1B1-guided statin selection and dosing |
| 6 | Chemotherapy Toxicity | Chemo Safety | DPYD/UGT1A1 screening for fluoropyrimidine/irinotecan toxicity |
| 7 | HLA Screening | HLA Screening | Pre-prescription HLA allele screening for drug hypersensitivity |
| 8 | Polypharmacy Review | Medication Review | Multi-drug phenoconversion detection and interaction analysis |

---

## 9. RAG Engine

**File:** `src/rag_engine.py` (799 lines)

### 9.1 PGxRAGEngine -- Core Methods

| Method | Purpose |
|--------|---------|
| `retrieve` | Embed query, parallel search 15 collections, expand, merge/rank, knowledge augment |
| `query` | Full RAG: retrieve + LLM synthesis (Claude Sonnet 4.6, 2048 tokens) |
| `query_stream` | Streaming variant: yields evidence dict then token-by-token LLM output |
| `find_related` | Cross-collection entity linking (e.g., "everything about CYP2D6") |
| `_embed_query` | BGE-small-en-v1.5 with retrieval prefix |
| `_search_all_collections` | ThreadPoolExecutor across 15 collections with weighted scores |
| `_merge_and_rank` | Deduplicate, sort by score descending, cap at 30 results |
| `_get_knowledge_context` | Extract structured context from pharmacogene/drug/HLA/phenoconversion knowledge |
| `_build_prompt` | Evidence by collection + knowledge context + citation instructions |

### 9.2 System Prompt

11 domains of expertise covering pharmacogene interpretation, drug-gene interactions, star allele nomenclature, diplotype-to-phenotype translation, CPIC/DPWG/FDA guidelines, HLA screening, phenoconversion, multi-gene interactions, population pharmacogenetics, dosing algorithms, and clinical implementation.

### 9.3 Collection Search Weights

Weights sum to 1.0 and are configurable via environment variables:

| Collection | Weight | Rationale |
|-----------|--------|-----------|
| Drug Guidelines | 0.14 | Highest: clinical recommendations are the primary output |
| Drug Interactions | 0.12 | Core PGx: gene-drug interaction evidence |
| Gene Reference | 0.10 | Allele definitions needed for every query |
| HLA Hypersensitivity | 0.10 | Safety-critical: never miss an HLA flag |
| Clinical Evidence | 0.08 | Published outcomes and validation studies |
| Phenoconversion | 0.08 | Drug-drug-gene interaction detection |
| Dosing Algorithms | 0.07 | Quantitative dose calculations |
| FDA Labels | 0.06 | Regulatory pharmacogenomic labeling |
| Population Data | 0.06 | Population-specific allele frequencies |
| Drug Alternatives | 0.05 | Therapeutic substitution options |
| Clinical Trials | 0.04 | Ongoing PGx research |
| Genomic Evidence | 0.03 | Shared variant data from pipeline |
| Patient Profiles | 0.03 | Example diplotype profiles |
| Implementation | 0.02 | Clinical program guidance |
| Education | 0.02 | Educational materials |

---

## 10. Agent Architecture

**File:** `src/agent.py` (588 lines)

### 10.1 Five-Phase Pipeline

```
Question -> Plan -> Search -> Evaluate -> Synthesize -> Report
```

| Phase | Method | Detail |
|-------|--------|--------|
| 1. Plan | `search_plan()` | Identify genes, drugs, phenotypes, workflow type |
| 2. Search | `rag.retrieve()` | Execute with workflow-based weight boosting |
| 3. Evaluate | `evaluate_evidence()` | Classify as sufficient/partial/insufficient |
| 4. Synthesize | `rag.query()` | LLM synthesis with clinical pipeline output |
| 5. Report | `generate_report()` | Structured PGx report with alerts and dosing |

### 10.2 Workflow Detection

The agent automatically detects 8 PGx workflow types from keyword analysis and routes to specialized processing pipelines.

---

## 11. Query Expansion

**File:** `src/query_expansion.py` (1,254 lines)

14 domain-specific expansion maps enrich queries with pharmacogenomic synonyms, brand/generic drug names, gene aliases, phenotype terms, and clinical concepts.

---

## 12. Data Models

**File:** `src/models.py` (616 lines)

Defines Pydantic models and enums for the PGx domain:

**Enums:** MetabolizerPhenotype, TransporterFunction, HLAStatus, EnzymeDeficiency, GuidelineBody, CPICLevel, ClinicalAction, AlertLevel, PGxWorkflowType

**Collection Models:** 14 Pydantic models mapping to Milvus collection schemas

**Query/Response Models:** PGxQuery, PGxResponse, AgentQuery, SearchHit, CrossCollectionResult, ComparativeResult, PGxAlert

---

## 13. Export System

**File:** `src/export.py` (1,307 lines)

Four export formats with PGx-specific enhancements:

| Format | Function | PGx Additions |
|--------|----------|---------------|
| Markdown | `export_markdown()` | Alert table, drug interaction matrix |
| JSON | `export_json()` | Pydantic serialization with alert objects |
| PDF | `export_pdf()` | Styled report with PGx Passport format |
| FHIR R4 | `export_fhir_r4()` | DiagnosticReport Bundle with LOINC 69548-6 PGx Observations |

---

## 14. API Reference

**Port:** 8107

### Core Endpoints (8)

| Method | Path | Description |
|--------|------|-------------|
| GET | `/` | Service root |
| GET | `/health` | Health check with collection/vector counts |
| GET | `/collections` | Collection names and record counts |
| POST | `/query` | Full RAG query (retrieve + LLM synthesis) |
| POST | `/search` | Evidence-only retrieval (no LLM) |
| POST | `/find-related` | Cross-collection entity linking |
| GET | `/knowledge/stats` | Knowledge graph statistics |
| GET | `/metrics` | Prometheus-compatible metrics |

### PGx Clinical Endpoints (7)

| Method | Path | Description |
|--------|------|-------------|
| POST | `/v1/pgx/drug-check` | Single drug PGx check with alerts |
| POST | `/v1/pgx/medication-review` | Polypharmacy medication review |
| POST | `/v1/pgx/warfarin-dose` | IWPC warfarin dosing calculation |
| POST | `/v1/pgx/hla-screen` | HLA hypersensitivity screening |
| POST | `/v1/pgx/phenoconversion` | Phenoconversion analysis |
| GET | `/v1/pgx/gene/{gene}` | Gene reference lookup |
| GET | `/v1/pgx/drug/{drug}` | Drug guideline lookup |

### Additional Endpoints

| Method | Path | Description |
|--------|------|-------------|
| POST | `/v1/events` | Event audit trail |
| POST | `/v1/reports` | Report generation |

---

## 15. UI Guide

**File:** `app/pgx_ui.py` (2,152 lines)
**Port:** 8507

### 10 Tabs

| # | Tab | Purpose |
|---|-----|---------|
| 1 | Dashboard | Overview with gene panel results, collection stats, knowledge graph |
| 2 | Drug Check | Single-drug PGx analysis with gene/phenotype/diplotype input |
| 3 | Medication Review | Multi-drug polypharmacy review with interaction detection |
| 4 | Warfarin Dosing | IWPC algorithm calculator with patient parameters |
| 5 | Chemo Safety | DPYD/UGT1A1 chemotherapy toxicity screening |
| 6 | HLA Screening | Pre-prescription HLA allele screening |
| 7 | Report Generator | Generate clinical PGx reports in MD/JSON/PDF/FHIR |
| 8 | Evidence Explorer | Browse and search across all 15 collections |
| 9 | Phenoconversion Modeler | Interactive CYP inhibitor/inducer phenoconversion analysis |
| 10 | Population Analytics | Population-specific allele frequency comparisons |

### UI Features

- NVIDIA dark theme (green #76B900 on dark background)
- Deep Research mode (autonomous agent pipeline)
- Conversation memory for follow-up queries (up to 3 prior exchanges)
- Collection-specific filtering in sidebar
- Citation relevance scoring (High >= 0.75, Medium >= 0.60, Low < 0.60)
- Export buttons for all 4 formats

---

## 16. Metrics & Monitoring

**File:** `src/metrics.py` (399 lines)

22 Prometheus metrics with `pgx_` prefix:

### Histograms (10)

| Metric | Labels | Description |
|--------|--------|-------------|
| `pgx_query_latency_seconds` | query_type | Query processing time |
| `pgx_evidence_count` | — | Evidence items per query |
| `pgx_cross_collection_query_latency_seconds` | query_type | Cross-collection query time |
| `pgx_cross_collection_results_count` | — | Cross-collection result count |
| `pgx_llm_api_latency_seconds` | provider, model | LLM API call latency |
| `pgx_embedding_latency_seconds` | — | Embedding generation time |
| `pgx_pipeline_stage_latency_seconds` | stage | Clinical pipeline stage time |
| `pgx_dosing_calculation_seconds` | algorithm | Dosing algorithm execution time |
| `pgx_hla_screening_seconds` | — | HLA screening execution time |
| `pgx_phenoconversion_check_seconds` | — | Phenoconversion check time |

### Counters (8)

| Metric | Labels | Description |
|--------|--------|-------------|
| `pgx_queries_total` | query_type | Total queries by type |
| `pgx_errors_total` | error_type | Total errors by type |
| `pgx_alerts_generated_total` | severity | Clinical alerts generated |
| `pgx_drug_checks_total` | — | Drug check requests |
| `pgx_hla_screens_total` | result | HLA screening requests |
| `pgx_dosing_calculations_total` | algorithm | Dosing calculations |
| `pgx_phenoconversion_detections_total` | — | Phenoconversion events detected |
| `pgx_export_total` | format | Export operations by format |

### Gauges (4)

| Metric | Description |
|--------|-------------|
| `pgx_collections_connected` | Number of connected collections |
| `pgx_total_vectors` | Total vectors across all collections |
| `pgx_last_ingest_timestamp` | Unix timestamp of last ingest run |
| `pgx_active_sessions` | Active UI sessions |

---

## 17. Scheduler

**File:** `src/scheduler.py` (232 lines)

APScheduler BackgroundScheduler with two recurring jobs:
- **PubMed PGx refresh**: Weekly (168 hours) fetch of new pharmacogenomics literature
- **ClinicalTrials.gov refresh**: Weekly fetch of PGx-related trial updates

Configurable via `PGX_INGEST_SCHEDULE_HOURS` and `PGX_INGEST_ENABLED` environment variables.

---

## 18. Configuration

**File:** `config/settings.py` (191 lines)

Pydantic BaseSettings with `PGX_` environment variable prefix. Key settings:

| Setting | Default | Description |
|---------|---------|-------------|
| `MILVUS_HOST` | localhost | Milvus server host |
| `MILVUS_PORT` | 19530 | Milvus server port |
| `LLM_MODEL` | claude-sonnet-4-6 | LLM model ID |
| `EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Embedding model |
| `EMBEDDING_DIMENSION` | 384 | Embedding vector dimension |
| `TOP_K_PER_COLLECTION` | 5 | Results per collection |
| `SCORE_THRESHOLD` | 0.4 | Minimum similarity score |
| `API_PORT` | 8107 | FastAPI server port |
| `STREAMLIT_PORT` | 8507 | Streamlit UI port |
| `INGEST_SCHEDULE_HOURS` | 168 | Ingest refresh interval |
| `MAX_CONVERSATION_CONTEXT` | 3 | Prior exchanges in context |

---

## 19. Docker Deployment

**File:** `docker-compose.yml`

6-service stack:

| Service | Image | Port | Purpose |
|---------|-------|------|---------|
| `milvus-etcd` | coreos/etcd:v3.5.5 | — | Milvus metadata store |
| `milvus-minio` | minio:RELEASE.2023-03-20 | — | Milvus object storage |
| `milvus-standalone` | milvus:v2.4-latest | 19530, 9091 | Vector database |
| `pgx-streamlit` | Custom | 8507 | Streamlit UI |
| `pgx-api` | Custom | 8107 | FastAPI REST server |
| `pgx-setup` | Custom | — | One-shot: create collections + seed data |

---

## 20. Testing

1,001 tests across 17 test files, all passing in 0.48 seconds.

| Test File | Coverage |
|-----------|----------|
| `test_knowledge.py` | 25 pharmacogenes, entity aliases, knowledge context functions |
| `test_collections.py` | Schema definitions, parallel search, weight normalization |
| `test_models.py` | Pydantic model validation, enum membership, serialization |
| `test_rag_engine.py` | Retrieval, synthesis, streaming, comparative analysis |
| `test_agent.py` | Plan generation, evidence evaluation, report formatting |
| `test_pgx_pipeline.py` | Star allele calling, phenotype translation, drug matching, alerts |
| `test_dosing.py` | IWPC warfarin, tacrolimus, fluoropyrimidine, thiopurine algorithms |
| `test_hla_screener.py` | All 15 HLA-drug associations, population risk, status assignment |
| `test_phenoconversion.py` | CYP inhibitor/inducer detection, phenotype shift modeling |
| `test_query_expansion.py` | All 14 expansion maps, deduplication, category separation |
| `test_export.py` | Markdown, JSON, PDF, FHIR R4 rendering with alerts |
| `test_metrics.py` | All 22 Prometheus metric registrations |
| `test_settings.py` | PGX_ prefix resolution, default values, type validation |
| `test_scheduler.py` | Job registration, interval configuration, no-op stub |

---

## 21. Port Map

| Port | Service | Protocol |
|------|---------|----------|
| 8507 | Streamlit UI | HTTP |
| 8107 | FastAPI REST API | HTTP |
| 19530 | Milvus gRPC | gRPC |
| 9091 | Milvus health/metrics | HTTP |

---

## 22. Tech Stack

| Layer | Technology |
|-------|-----------|
| Compute | NVIDIA DGX Spark (GB10 GPU, 128GB unified memory) |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Embeddings | BGE-small-en-v1.5 (384-dim, HuggingFace) |
| Vector DB | Milvus 2.4 (IVF_FLAT, COSINE) |
| Backend | FastAPI, Pydantic, uvicorn |
| Frontend | Streamlit |
| Clinical Logic | Custom Python (star allele calling, dosing algorithms, HLA screening) |
| Export | reportlab (PDF), FHIR R4 (JSON) |
| Monitoring | Prometheus (prometheus_client), Grafana |
| Scheduling | APScheduler |
| Container | Docker, Docker Compose |
| Testing | pytest |

---

## 23. Future Roadmap

1. **CYP2D6 structural variant calling** -- Integration with Stargazer/Cyrius for copy number and hybrid allele detection
2. **EHR integration** -- HL7 FHIR CDS Hooks for real-time clinical decision support in Epic/Cerner
3. **Multi-language support** -- Spanish, Mandarin, and Arabic for patient-facing PGx reports
4. **Federated learning** -- Privacy-preserving model updates across institutional PGx databases
5. **Expanded dosing algorithms** -- CYP2C19-guided voriconazole, SLCO1B1-guided simvastatin
6. **Pharmacovigilance integration** -- FDA FAERS adverse event signal detection
7. **Pre-emptive panel workflows** -- End-to-end integration with genomics pipeline VCF output
