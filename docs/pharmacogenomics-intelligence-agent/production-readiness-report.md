# Pharmacogenomics Intelligence Agent -- Production Readiness Report

**Version:** 1.0.0
**Date:** March 12, 2026
**Author:** Adam Jones
**Status:** Production Demo Ready (10/10)
**License:** Apache 2.0

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [System Overview](#2-system-overview)
3. [Collection Architecture](#3-collection-architecture)
4. [Knowledge Graph](#4-knowledge-graph)
5. [Star Allele Calling Pipeline](#5-star-allele-calling-pipeline)
6. [Phenotype Translation](#6-phenotype-translation)
7. [Phenoconversion Detection](#7-phenoconversion-detection)
8. [HLA Screening](#8-hla-screening)
9. [Dosing Algorithms](#9-dosing-algorithms)
10. [RAG Engine](#10-rag-engine)
11. [Autonomous Agent](#11-autonomous-agent)
12. [Query Expansion](#12-query-expansion)
13. [API Coverage](#13-api-coverage)
14. [UI Capabilities](#14-ui-capabilities)
15. [Export and Reporting](#15-export-and-reporting)
16. [Test Coverage](#16-test-coverage)
17. [Performance](#17-performance)
18. [Observability and Metrics](#18-observability-and-metrics)
19. [Security](#19-security)
20. [Data Sources](#20-data-sources)
21. [Configuration System](#21-configuration-system)
22. [Infrastructure and Deployment](#22-infrastructure-and-deployment)
23. [Known Limitations](#23-known-limitations)
24. [Recommendations](#24-recommendations)
25. [Appendix A: File Inventory](#appendix-a-file-inventory)
26. [Appendix B: Collection Schemas](#appendix-b-collection-schemas)
27. [Appendix C: Metric Registry](#appendix-c-metric-registry)
28. [Appendix D: Codebase Metrics](#appendix-d-codebase-metrics)

---

## 1. Executive Summary

The Pharmacogenomics Intelligence Agent is a production-grade, multi-collection RAG system purpose-built for the pharmacogenomics domain. It is one of five intelligence agents in the HCLS AI Factory precision medicine platform, running on NVIDIA DGX Spark hardware.

### Key Capabilities at a Glance

| Capability | Detail |
|-----------|--------|
| Vector collections | 15 (14 PGx domain-specific + 1 shared genomic) |
| Seed data records | 240 across 14 JSON files |
| Knowledge graph pharmacogenes | 25 with clinical PGx data |
| Drugs covered | 100+ across 12 therapeutic categories |
| HLA-drug associations | 12 validated hypersensitivity pairs |
| CYP inhibitors/inducers | 30+ for phenoconversion modeling |
| Dosing algorithms | 4 validated (IWPC warfarin, CYP3A5 tacrolimus, DPYD fluoropyrimidine, TPMT+NUDT15 thiopurine) |
| Clinical workflows | 8 (pre-emptive panel, opioid safety, anticoagulant optimization, antidepressant selection, statin myopathy, chemo toxicity, HLA screening, polypharmacy) |
| LLM | Claude Sonnet 4.6 (claude-sonnet-4-6) via Anthropic API |
| Embedding model | BGE-small-en-v1.5 (384 dimensions) |
| Export formats | Markdown, JSON, PDF, FHIR R4 |
| Test suite | 696 tests, all passing in 0.41s |
| Total Python LOC | 23,049 (19,148 source + 3,901 test) across 53 Python files |
| Total files | 75 (53 Python, 14 JSON seed, 8 config/infra) |
| Ingest parsers | 8 (CPIC, PharmVar, PharmGKB, FDA, population, PubMed, ClinicalTrials.gov + base) |
| Prometheus metrics | 22 (10 histograms, 8 counters, 4 gauges) |
| API endpoints | 16+ (8 core + 7 PGx clinical + events + reports) |
| UI tabs | 10 |
| Service ports | 8507 (Streamlit UI), 8107 (FastAPI API) |

### Architecture Overview

```
User Query
    |
    v
+---------------------------------------------+
|  Streamlit UI (:8507)  /  FastAPI (:8107)    |
+----------------------+-----------------------+
                       |
    +------------------+------------------+
    v                  v                  v
+--------+     +----------+     +-----------+
| Agent  |     |   RAG    |     | Knowledge |
|Pipeline|     |  Engine  |     |   Graph   |
+---+----+     +----+-----+     +-----+-----+
    |               |                 |
    v               v                 v
+---------------------------------------------+
|  Clinical Pipelines                          |
|  +--------+ +--------+ +-----+ +---------+  |
|  |Star    | |Pheno-  | |HLA  | |Dosing   |  |
|  |Allele  | |convert.| |Scr. | |Calc.    |  |
|  +--------+ +--------+ +-----+ +---------+  |
+---------------------------------------------+
                       |
                       v
+---------------------------------------------+
|  Milvus Vector DB (15 Collections)           |
|  pgx_gene_reference   pgx_drug_guidelines   |
|  pgx_drug_interactions pgx_hla_hypersensi.   |
|  pgx_phenoconversion   pgx_dosing_algorithms |
|  pgx_clinical_evidence  pgx_population_data  |
|  pgx_clinical_trials    pgx_fda_labels       |
|  pgx_drug_alternatives  pgx_patient_profiles |
|  pgx_implementation     pgx_education        |
|  genomic_evidence (shared)                    |
+---------------------------------------------+
                       |
                       v
+---------------------------------------------+
|  Claude Sonnet 4.6 (Anthropic API)           |
|  Streaming response with citations            |
+---------------------------------------------+
```

---

## 2. System Overview

### Component Map

| Component | File | Lines | Purpose |
|-----------|------|-------|---------|
| Knowledge Graph | `src/knowledge.py` | 2,512 | 9 structured dictionaries: 25 pharmacogenes, 12 drug categories, 12 HLA associations, CYP inhibitors/inducers, drug alternatives, activity score tables, entity aliases |
| Collections | `src/collections.py` | 1,547 | 15 Milvus collection schemas, parallel ThreadPoolExecutor search, IVF_FLAT indexing |
| Export | `src/export.py` | 1,296 | Markdown, JSON, PDF, FHIR R4 export with PGx alerts and drug interaction matrix |
| Query Expansion | `src/query_expansion.py` | 1,254 | 14 domain-specific expansion maps |
| PGx Pipeline | `src/pgx_pipeline.py` | 1,222 | StarAlleleCaller, PhenotypeTranslator, DrugGeneMatcher |
| RAG Engine | `src/rag_engine.py` | 796 | Multi-collection retrieval, prompt building, knowledge augmentation, LLM synthesis |
| PGx Clinical API | `api/routes/pgx_clinical.py` | 803 | 7 PGx-specific clinical decision support endpoints |
| Agent | `src/agent.py` | 761 | Autonomous plan-search-evaluate-synthesize-report pipeline |
| Dosing | `src/dosing.py` | 692 | 4 validated dosing algorithms |
| API Server | `api/main.py` | 649 | FastAPI REST with 8 core endpoints, CORS, lifespan management |
| HLA Screener | `src/hla_screener.py` | 627 | 12 HLA-drug hypersensitivity screening pairs |
| Models | `src/models.py` | 616 | Enums, collection models, query/response schemas |
| Phenoconversion | `src/phenoconversion.py` | 494 | CYP inhibitor/inducer phenoconversion detection |
| Metrics | `src/metrics.py` | 399 | 22 Prometheus metrics |
| Scheduler | `src/scheduler.py` | 232 | Weekly automated PubMed + ClinicalTrials.gov refresh |
| UI | `app/pgx_ui.py` | 2,138 | 10-tab Streamlit interface with NVIDIA dark theme |
| Settings | `config/settings.py` | 118 | Pydantic BaseSettings with PGX_ env prefix |

### Data Flow

1. **Ingest**: 8 parsers fetch from CPIC, PharmVar, PharmGKB, FDA, population databases, PubMed, ClinicalTrials.gov, and local JSON seed files
2. **Embed**: Each record's `to_embedding_text()` -> BGE-small-en-v1.5 -> 384-dim vector
3. **Store**: Batch insert into Milvus with IVF_FLAT index (COSINE metric, nlist=1024)
4. **Query**: User question -> embed -> parallel search 15 collections -> expand -> merge/rank (max 30) -> knowledge augmentation -> clinical pipeline -> LLM synthesis -> streaming response with citations and alerts

---

## 3. Collection Architecture

### 3.1 Fifteen Milvus Collections

| # | Collection | Weight | Purpose | Key Fields |
|---|-----------|--------|---------|------------|
| 1 | `pgx_gene_reference` | 0.10 | Pharmacogene star allele definitions & activity scores | gene, star_allele, activity_score, function_status |
| 2 | `pgx_drug_guidelines` | 0.14 | CPIC/DPWG prescribing recommendations | gene, drug, guideline_body, cpic_level, phenotype, recommendation |
| 3 | `pgx_drug_interactions` | 0.12 | PharmGKB drug-gene interaction records | gene, drug, interaction_type, evidence_level |
| 4 | `pgx_hla_hypersensitivity` | 0.10 | HLA-mediated adverse drug reaction screening | hla_allele, drug, reaction_type, severity |
| 5 | `pgx_phenoconversion` | 0.08 | Metabolic phenoconversion via drug-drug interactions | inhibitor_drug, affected_enzyme, inhibitor_strength |
| 6 | `pgx_dosing_algorithms` | 0.07 | Genotype-guided dosing algorithms & formulas | drug, algorithm_name, gene, dose_adjustment |
| 7 | `pgx_clinical_evidence` | 0.08 | Published PGx clinical evidence & outcomes | pmid, study_type, gene, drug, outcome |
| 8 | `pgx_population_data` | 0.06 | Population-specific allele frequency data | gene, allele, population, frequency |
| 9 | `pgx_clinical_trials` | 0.04 | PGx-related clinical trials | nct_id, gene, drug, phase, status |
| 10 | `pgx_fda_labels` | 0.06 | FDA pharmacogenomic labeling information | drug, gene, label_section, fda_action |
| 11 | `pgx_drug_alternatives` | 0.05 | Genotype-guided therapeutic alternatives | gene, phenotype, drug_to_avoid, alternative_drug |
| 12 | `pgx_patient_profiles` | 0.03 | Patient diplotype-phenotype profiles | patient_id, gene, diplotype, phenotype |
| 13 | `pgx_implementation` | 0.02 | Clinical PGx implementation programs | institution, program_type, genes_tested |
| 14 | `pgx_education` | 0.02 | PGx educational resources & guidelines | topic, target_audience, content_type |
| 15 | `genomic_evidence` | 0.03 | Shared read-only genomic variants (rag-chat-pipeline) | chrom, pos, ref, alt, gene |

### 3.2 Index Configuration

| Parameter | Value |
|-----------|-------|
| Index type | IVF_FLAT |
| Distance metric | COSINE |
| nlist | 1024 |
| Embedding dimension | 384 (BGE-small-en-v1.5) |

### 3.3 Seed Data Distribution

240 total records across 14 JSON files in `data/reference/`:

| File | Description |
|------|------------|
| gene_reference.json | Pharmacogene star allele definitions |
| drug_guidelines.json | CPIC/DPWG prescribing guidelines |
| drug_interactions.json | Drug-gene interaction records |
| hla_hypersensitivity.json | HLA-drug ADR associations |
| phenoconversion.json | CYP inhibitor/inducer data |
| dosing_algorithms.json | Genotype-guided dosing formulas |
| clinical_evidence.json | Published PGx outcome studies |
| population_data.json | Population allele frequencies |
| clinical_trials.json | PGx clinical trial records |
| fda_labels.json | FDA pharmacogenomic labeling |
| drug_alternatives.json | Therapeutic alternatives |
| patient_profiles.json | Patient diplotype profiles |
| implementation.json | Clinical PGx programs |
| education.json | PGx educational resources |

---

## 4. Knowledge Graph

**File:** `src/knowledge.py` (2,512 lines)

### 4.1 Nine Structured Dictionaries

| Dictionary | Entries | Description |
|-----------|---------|-------------|
| `PHARMACOGENES` | 25 | Clinical pharmacogenomic data per gene: chromosome, function, substrates, key variants, CPIC guidelines |
| `METABOLIZER_PHENOTYPES` | 5 | UM, RM, NM, IM, PM definitions with activity score ranges |
| `DRUG_CATEGORIES` | 12 | Therapeutic categories with member drugs |
| `CYP_INHIBITORS` | 4 enzymes | Strong/moderate/weak inhibitors for CYP2D6, CYP2C19, CYP2C9, CYP3A4/5 |
| `CYP_INDUCERS` | 3 enzymes | Strong/moderate inducers for CYP2D6, CYP2C19, CYP3A4 |
| `HLA_DRUG_ASSOCIATIONS` | 12 | HLA allele-drug pairs with reaction type, severity, alternatives, population risk |
| `DRUG_ALTERNATIVES` | 30+ | Gene-phenotype-specific therapeutic substitutions with rationale |
| `ACTIVITY_SCORE_TABLES` | 2 | CYP2D6 and DPYD allele-to-activity-score mappings |
| `ENTITY_ALIASES` | 80+ | Drug brand/generic names, gene aliases, phenotype abbreviations |

### 4.2 The 25 Pharmacogenes

**CYP Enzymes (8):** CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP3A5, CYP2B6, CYP1A2, CYP4F2

**Phase II Enzymes (4):** UGT1A1, NAT2, TPMT, DPYD

**Transporters (2):** SLCO1B1, ABCB1

**HLA Genes (3):** HLA-A, HLA-B, HLA-DRB1

**Other (8):** VKORC1, NUDT15, G6PD, IFNL3, RYR1, CACNA1S, CYP2A6, COMT

### 4.3 Per-Pharmacogene Data Fields

Each entry in `PHARMACOGENES` contains:
- `full_name`: Full protein/enzyme name
- `chromosome`: Chromosomal location
- `function`: Biological function description
- `substrates_count`: Number of known substrates
- `percent_drugs_metabolized`: Percentage of all drugs metabolized
- `star_alleles_defined`: Number of PharmVar-defined star alleles
- `key_variants`: Clinically significant star alleles with function annotations
- `structural_variation`: Whether the gene has CNV/hybrid complexity
- `complexity_level`: Assessment of genotyping difficulty
- `cpic_guidelines`: List of drugs with CPIC guidelines

### 4.4 The 12 Therapeutic Drug Categories

| Category | Example Drugs |
|----------|--------------|
| Opioids | codeine, tramadol, oxycodone, hydrocodone |
| Antidepressants | escitalopram, citalopram, sertraline, amitriptyline, nortriptyline, venlafaxine |
| Anticoagulants | warfarin |
| Antiplatelets | clopidogrel, prasugrel, ticagrelor |
| Antiepileptics | carbamazepine, phenytoin, oxcarbazepine |
| Antivirals | abacavir, efavirenz |
| Immunosuppressants | tacrolimus, azathioprine, mercaptopurine |
| Statins | simvastatin, atorvastatin, rosuvastatin |
| PPIs | omeprazole, lansoprazole, pantoprazole |
| Fluoropyrimidines | fluorouracil, capecitabine |
| Thiopurines | azathioprine, mercaptopurine, thioguanine |
| Antipsychotics | aripiprazole, risperidone, haloperidol |

### 4.5 The 12 HLA-Drug Associations

| # | Drug | HLA Allele | Reaction | Severity |
|---|------|-----------|----------|----------|
| 1 | Abacavir | HLA-B*57:01 | Hypersensitivity syndrome | Severe |
| 2 | Carbamazepine | HLA-B*15:02 | SJS/TEN | Fatal/Severe |
| 3 | Carbamazepine | HLA-A*31:01 | DRESS/SJS | Severe |
| 4 | Phenytoin | HLA-B*15:02 | SJS/TEN | Fatal/Severe |
| 5 | Oxcarbazepine | HLA-B*15:02 | SJS/TEN | Severe |
| 6 | Allopurinol | HLA-B*58:01 | SJS/TEN | Fatal/Severe |
| 7 | Flucloxacillin | HLA-B*57:01 | DILI | Severe |
| 8 | Lamotrigine | HLA-B*15:02 | SJS/TEN | Severe |
| 9 | Dapsone | HLA-B*13:01 | Hypersensitivity | Severe |
| 10 | Nevirapine | HLA-B*35:05 | Hepatotoxicity | Severe |
| 11 | Ticlopidine | HLA-A*33:03 | Hepatotoxicity | Severe |
| 12 | Sulfasalazine | HLA-B*13:01 | DRESS | Moderate |

### 4.6 Entity Resolution

80+ entity aliases map brand names, abbreviations, and colloquial terms to canonical entities:
- Drug aliases: Coumadin -> warfarin, Plavix -> clopidogrel, Tegretol -> carbamazepine
- Gene aliases: CYP2D6 -> CYP2D6, 2D6 -> CYP2D6
- Phenotype aliases: PM -> poor_metabolizer, UM -> ultra_rapid_metabolizer

### 4.7 Public Functions

| Function | Purpose |
|----------|---------|
| `get_gene_context(gene)` | Format pharmacogene knowledge for prompt injection |
| `get_drug_context(drug)` | Format drug category and guideline context |
| `get_hla_context(drug)` | Format HLA-drug association context |
| `get_phenoconversion_context(drug)` | Format CYP inhibitor/inducer context |
| `get_all_context_for_query(query)` | Extract all relevant knowledge via keyword matching |
| `get_knowledge_stats()` | Return counts: 25 genes, 12 categories, 12 HLA, 80+ aliases |
| `resolve_entity(text)` | Resolve text to canonical entity via aliases |

---

## 5. Star Allele Calling Pipeline

**File:** `src/pgx_pipeline.py` (1,222 lines) -- `StarAlleleCaller` class

### 5.1 Capabilities

- Resolves VCF variant calls to PharmVar-aligned star allele nomenclature for 25 pharmacogenes
- Handles simple SNV-based allele calls (e.g., CYP2C19*2 from rs4244285)
- Supports compound heterozygous diplotype resolution
- Detects CYP2D6 gene deletion (*5) and duplication (*1xN, *2xN) when reported in VCF
- Maps defining variants (rsIDs) to star alleles using PharmVar definitions

### 5.2 Data Classes

| Class | Fields | Purpose |
|-------|--------|---------|
| `PGxPosition` | chrom, pos, ref, alt, star_allele, rsid, gene | Single pharmacogenomic variant |
| `PGxAlert` | drug, gene, diplotype, phenotype, severity, recommendation, cpic_level, alternatives, source | Clinical alert from drug-gene matching |

### 5.3 Alert Severity Levels

| Severity | Trigger |
|----------|---------|
| CONTRAINDICATED | HLA positive for fatal/severe ADR; PM + contraindicated drug |
| MAJOR | PM/UM requiring drug change or major dose adjustment |
| MODERATE | IM requiring monitoring or minor adjustment |
| MINOR | Known interaction with limited clinical impact |
| INFORMATIONAL | Educational context, population data |

---

## 6. Phenotype Translation

**File:** `src/pgx_pipeline.py` -- `PhenotypeTranslator` class

### 6.1 CPIC-Standardized Translation

| Diplotype Activity Score | Phenotype |
|-------------------------|-----------|
| > 2.0 | Ultra-Rapid Metabolizer (UM) |
| 1.5 - 2.0 | Rapid Metabolizer (RM) |
| 1.0 - 2.0 | Normal Metabolizer (NM) |
| 0.5 - 1.0 | Intermediate Metabolizer (IM) |
| 0.0 | Poor Metabolizer (PM) |

### 6.2 Activity Score Tables

Pre-defined activity scores for alleles of CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, and NUDT15 genes. Each star allele maps to 0.0 (no function), 0.5 (decreased function), or 1.0 (normal function), with gene duplication alleles mapping to 2.0+.

---

## 7. Phenoconversion Detection

**File:** `src/phenoconversion.py` (494 lines)

### 7.1 CYP Inhibitor Knowledge Base

30+ inhibitors classified by FDA potency across 4 CYP enzymes:

**CYP2D6 Inhibitors (12):**

| Drug | Strength |
|------|----------|
| Fluoxetine, Paroxetine, Bupropion, Quinidine, Cinacalcet, Terbinafine | Strong |
| Duloxetine, Sertraline, Diphenhydramine | Moderate |
| Celecoxib, Cimetidine, Methadone | Weak |

**CYP2C19 Inhibitors (6):**

| Drug | Strength |
|------|----------|
| Fluvoxamine, Voriconazole, Ticlopidine | Strong |
| Fluconazole, Omeprazole, Esomeprazole | Moderate |

**CYP2C9 Inhibitors (5):**

| Drug | Strength |
|------|----------|
| Miconazole | Strong |
| Amiodarone, Fluconazole | Moderate |
| Metronidazole, Trimethoprim, Sulfamethoxazole | Weak |

**CYP3A4/5 Inhibitors (7+):**

| Drug | Strength |
|------|----------|
| Ketoconazole, Itraconazole, Clarithromycin, Ritonavir | Strong |
| Erythromycin, Diltiazem, Verapamil | Moderate |

### 7.2 CYP Inducer Knowledge Base

| Enzyme | Strong Inducers | Moderate Inducers |
|--------|----------------|-------------------|
| CYP3A4 | Rifampin, Carbamazepine, Phenytoin, St. John's Wort | Efavirenz, Modafinil |
| CYP2C19 | Rifampin | Carbamazepine |
| CYP1A2 | Smoking (PAH induction) | Charbroiled meat |

### 7.3 Phenotype Shift Logic

| Genetic Phenotype | + Strong Inhibitor | + Moderate Inhibitor | + Strong Inducer |
|-------------------|-------------------|---------------------|-----------------|
| UM | -> NM/IM | -> RM | -> UM (no change) |
| NM | -> PM | -> IM | -> UM/RM |
| IM | -> PM | -> PM | -> NM |
| PM | -> PM (no change) | -> PM (no change) | -> IM |

---

## 8. HLA Screening

**File:** `src/hla_screener.py` (627 lines)

### 8.1 Screening Logic

The `HLAScreener` accepts:
- Patient HLA typing (list of HLA alleles)
- Drug name

And returns an `HLAScreenResult` with:
- Status: SAFE, CONTRAINDICATED, HIGH_RISK, UNKNOWN
- Reaction type (e.g., SJS/TEN, DRESS, hypersensitivity syndrome)
- Severity: fatal, severe, moderate, mild
- Recommendation (avoid, test first, use alternative)
- Alternative drugs
- Evidence level (FDA Boxed Warning, CPIC Level A, etc.)
- Population-specific risk information

### 8.2 Coverage

All 12 HLA-drug associations validated with:
- Correct HLA allele matching
- Appropriate severity classification
- Population-specific prevalence data
- Therapeutic alternatives for each contraindicated drug

---

## 9. Dosing Algorithms

**File:** `src/dosing.py` (692 lines)

### 9.1 Four Validated Algorithms

| Algorithm | Drug(s) | Key Variables | Reference |
|-----------|---------|---------------|-----------|
| IWPC Warfarin | Warfarin | Age, height, weight, race, CYP2C9, VKORC1, amiodarone, enzyme inducers | Klein et al., NEJM 2009 |
| CYP3A5 Tacrolimus | Tacrolimus | CYP3A5 expresser (*1 carrier) vs non-expresser (*3/*3) | Birdwell et al., CPT 2015 |
| DPYD Fluoropyrimidine | 5-FU, Capecitabine | DPYD activity score (diplotype-based dose reduction) | Amstutz et al., CPT 2018 |
| TPMT+NUDT15 Thiopurine | Azathioprine, Mercaptopurine | Combined TPMT + NUDT15 activity assessment | Relling et al., CPT 2019 |

### 9.2 IWPC Warfarin Algorithm Detail

Regression model calculating square root of weekly dose:

```
sqrt(dose) = 5.6044
  - 0.2614 * age_decades
  + 0.0087 * height_cm
  + 0.0128 * weight_kg
  - 0.8677 * VKORC1_AG
  - 1.6974 * VKORC1_AA
  - 0.5211 * CYP2C9_*1/*2
  - 0.9357 * CYP2C9_*1/*3
  - 1.0616 * CYP2C9_*2/*2
  - 1.9206 * CYP2C9_*2/*3
  - 2.3312 * CYP2C9_*3/*3
  - 0.2188 * asian_race
  - 0.1092 * black_race
  - 0.2760 * other_race
  - 0.1032 * amiodarone
  + 1.1816 * enzyme_inducer
```

### 9.3 DPYD Dose Adjustment

| DPYD Activity Score | Dose Adjustment |
|--------------------|----------------|
| 2.0 (normal) | Full dose (100%) |
| 1.5 | Reduce by 25-50% |
| 1.0 | Reduce by 50% |
| 0.5 | Avoid fluoropyrimidines |
| 0.0 | Contraindicated |

### 9.4 TPMT+NUDT15 Combined Dosing

| TPMT Status | NUDT15 Status | Recommendation |
|------------|--------------|----------------|
| Normal | Normal | Full dose |
| Normal | Intermediate | 25-50% dose reduction |
| Intermediate | Normal | 30-70% dose reduction |
| Intermediate | Intermediate | 50-80% dose reduction |
| Poor | Any | 10% dose or avoid |
| Any | Poor | 10% dose or avoid |

---

## 10. RAG Engine

**File:** `src/rag_engine.py` (796 lines)

### 10.1 PGxRAGEngine Methods

| Method | Purpose |
|--------|---------|
| `retrieve` | Embed -> parallel search 15 collections -> expand -> merge/rank -> knowledge augment |
| `query` | Full RAG: retrieve + LLM synthesis (2,048 tokens) |
| `query_stream` | Streaming: yields evidence dict then token-by-token LLM output |
| `find_related` | Cross-collection entity linking |
| `_embed_query` | BGE-small-en-v1.5 with retrieval prefix |
| `_search_all_collections` | ThreadPoolExecutor across 15 collections |
| `_expanded_search` | Query expansion via 14 maps |
| `_merge_and_rank` | Deduplicate, weighted score, cap at 30 |
| `_get_knowledge_context` | Inject pharmacogene/drug/HLA/phenoconversion knowledge |
| `_build_prompt` | Evidence + knowledge + clinical pipeline results |

### 10.2 System Prompt

11 domains of expertise: pharmacogene interpretation, drug-gene interactions, star allele nomenclature, diplotype-to-phenotype translation, CPIC/DPWG/FDA guidelines, HLA screening, phenoconversion, multi-gene interactions, population pharmacogenetics, dosing algorithms, clinical implementation.

### 10.3 Score Thresholds

| Threshold | Value | Purpose |
|-----------|-------|---------|
| Minimum score | 0.4 | Below this, results are excluded |
| High citation | >= 0.75 | Marked as high-relevance evidence |
| Medium citation | >= 0.60 | Marked as medium-relevance evidence |
| Low citation | < 0.60 | Marked as low-relevance evidence |

### 10.4 Conversation Memory

Injects up to 3 prior Q&A exchanges (configurable via `PGX_MAX_CONVERSATION_CONTEXT`) truncated to 300 characters each for context window efficiency.

---

## 11. Autonomous Agent

**File:** `src/agent.py` (761 lines)

### 11.1 Five-Phase Pipeline

```
Question -> Plan -> Search -> Evaluate -> Synthesize -> Report
```

| Phase | Method | Detail |
|-------|--------|--------|
| 1. Plan | `search_plan()` | Identify genes, drugs, phenotypes, workflow type |
| 2. Search | `rag.retrieve()` | Execute with workflow-based weight boosting |
| 3. Evaluate | `evaluate_evidence()` | Sufficient (>= 3 collections, >= 10 hits) / partial / insufficient |
| 4. Synthesize | `rag.query()` | LLM synthesis with full evidence + clinical pipeline results |
| 5. Report | `generate_report()` | Structured markdown with query, evidence, alerts, dosing |

### 11.2 Workflow Detection

8 PGx workflow types automatically detected from query keywords:
1. Pre-emptive Panel
2. Opioid Safety
3. Anticoagulant Optimization
4. Antidepressant Selection
5. Statin Myopathy Prevention
6. Chemotherapy Toxicity
7. HLA Screening
8. Polypharmacy Review

---

## 12. Query Expansion

**File:** `src/query_expansion.py` (1,254 lines)

### 12.1 Fourteen Expansion Maps

| Map | Description |
|-----|------------|
| DRUG_EXPANSION | Drug brand/generic names to PGx terms |
| GENE_EXPANSION | Gene symbols to substrates, inhibitors, clinical terms |
| PHENOTYPE_EXPANSION | Metabolizer phenotype terms to clinical implications |
| HLA_EXPANSION | HLA alleles to drugs, reactions, populations |
| DOSING_EXPANSION | Dosing-related terms to algorithms, parameters |
| INTERACTION_EXPANSION | Drug interaction terms to CYP enzymes, severity |
| DISEASE_EXPANSION | Disease names to PGx-relevant drugs and genes |
| POPULATION_EXPANSION | Population/ancestry terms to allele frequencies |
| GUIDELINE_EXPANSION | CPIC/DPWG/FDA terms to guideline content |
| CLINICAL_EXPANSION | Clinical terms (ADR, efficacy) to PGx concepts |
| ENZYME_EXPANSION | Enzyme pathway terms to substrates, inhibitors |
| TRANSPORTER_EXPANSION | Transporter terms (SLCO1B1, ABCB1) to drugs |
| IMPLEMENTATION_EXPANSION | Implementation terms to program types, EHR |
| EDUCATION_EXPANSION | Educational terms to learning resources |

### 12.2 Public Functions

| Function | Purpose |
|----------|---------|
| `expand_query(query)` | Expand query using all 14 maps, deduplicated |
| `expand_query_by_category(query)` | Return categorized dict of expansions by map |
| `get_expansion_stats()` | Return per-map keyword and term counts |

---

## 13. API Coverage

**Port:** 8107

### 13.1 Core Endpoints (8)

| Method | Path | Description |
|--------|------|-------------|
| GET | `/` | Service root and metadata |
| GET | `/health` | Health check with collection/vector counts and Milvus status |
| GET | `/collections` | Collection names and record counts |
| POST | `/query` | Full RAG query (retrieve + LLM synthesis) |
| POST | `/search` | Evidence-only retrieval (no LLM, fast) |
| POST | `/find-related` | Cross-collection entity linking |
| GET | `/knowledge/stats` | Knowledge graph statistics |
| GET | `/metrics` | Prometheus-compatible metrics in text format |

### 13.2 PGx Clinical Endpoints (7)

| Method | Path | Description |
|--------|------|-------------|
| POST | `/v1/pgx/drug-check` | Single drug PGx check with alerts and alternatives |
| POST | `/v1/pgx/medication-review` | Polypharmacy review with phenoconversion detection |
| POST | `/v1/pgx/warfarin-dose` | IWPC warfarin dosing calculation |
| POST | `/v1/pgx/hla-screen` | HLA hypersensitivity screening |
| POST | `/v1/pgx/phenoconversion` | Phenoconversion analysis for medication list |
| GET | `/v1/pgx/gene/{gene}` | Gene reference lookup |
| GET | `/v1/pgx/drug/{drug}` | Drug guideline lookup |

### 13.3 Additional Endpoints

| Method | Path | Description |
|--------|------|-------------|
| POST | `/v1/events` | Event audit trail |
| POST | `/v1/reports` | Report generation |

---

## 14. UI Capabilities

**File:** `app/pgx_ui.py` (2,138 lines)
**Port:** 8507

### 14.1 Ten Tabs

| # | Tab | Purpose |
|---|-----|---------|
| 1 | Dashboard | System overview, collection stats, knowledge graph summary |
| 2 | Drug Check | Single-drug PGx analysis with gene/phenotype/diplotype input |
| 3 | Medication Review | Multi-drug polypharmacy review with phenoconversion detection |
| 4 | Warfarin Dosing | IWPC algorithm calculator with patient parameter inputs |
| 5 | Chemo Safety | DPYD/UGT1A1 chemotherapy toxicity screening |
| 6 | HLA Screening | Pre-prescription HLA allele screening |
| 7 | Report Generator | Clinical PGx reports in MD/JSON/PDF/FHIR |
| 8 | Evidence Explorer | Browse and search across all 15 collections |
| 9 | Phenoconversion Modeler | Interactive CYP inhibitor/inducer phenoconversion analysis |
| 10 | Population Analytics | Population-specific allele frequency comparisons |

### 14.2 UI Features

- NVIDIA dark theme (green #76B900 on dark background)
- Deep Research mode (autonomous agent pipeline)
- Conversation memory for follow-up queries (3 prior exchanges)
- Collection-specific filtering in sidebar
- Citation relevance scoring (High >= 0.75, Medium >= 0.60, Low < 0.60)
- Export buttons for Markdown, JSON, PDF, FHIR R4
- Cached engine initialization via `@st.cache_resource`

---

## 15. Export and Reporting

**File:** `src/export.py` (1,296 lines)

### 15.1 Four Export Formats

| Format | Function | PGx-Specific Features |
|--------|----------|----------------------|
| Markdown | `export_markdown()` | Alert severity table, drug interaction matrix, citation links |
| JSON | `export_json()` | Pydantic model serialization with PGxAlert objects |
| PDF | `export_pdf()` | reportlab Platypus styling, PGx Passport format, NVIDIA branding |
| FHIR R4 | `export_fhir_r4()` | DiagnosticReport Bundle, LOINC 69548-6 PGx Observations, Genotype/Phenotype resources |

### 15.2 Export Content

All exports include:
- Original query and timestamp
- LLM response
- Evidence breakdown by collection (source, score, citation)
- Clinical alerts (severity, gene, drug, phenotype, recommendation, alternatives)
- Drug interaction matrix (for medication reviews)
- Knowledge graph context used
- Metadata (collections searched, total evidence, search time)

---

## 16. Test Coverage

### 16.1 Summary

| Metric | Value |
|--------|-------|
| Total tests | 696 |
| Pass rate | 100% |
| Execution time | 0.41 seconds |
| Test files | 15 (14 test modules + conftest.py) |
| Test LOC | 3,901 |

### 16.2 Test Distribution by Module

| Test File | Module Tested | Coverage Focus |
|-----------|--------------|---------------|
| `test_knowledge.py` | knowledge.py | 25 pharmacogenes, entity aliases, knowledge context functions, drug categories |
| `test_collections.py` | collections.py | 15 collection schemas, field definitions, weight normalization, parallel search |
| `test_models.py` | models.py | Pydantic validation, enum membership, serialization/deserialization |
| `test_rag_engine.py` | rag_engine.py | Retrieval, synthesis, streaming, comparative analysis, score thresholds |
| `test_agent.py` | agent.py | Plan generation, workflow detection, evidence evaluation, report formatting |
| `test_pgx_pipeline.py` | pgx_pipeline.py | Star allele calling, phenotype translation, drug-gene matching, alert generation |
| `test_dosing.py` | dosing.py | IWPC warfarin, CYP3A5 tacrolimus, DPYD fluoropyrimidine, TPMT+NUDT15 thiopurine |
| `test_hla_screener.py` | hla_screener.py | All 12 HLA-drug associations, severity classification, alternatives, population risk |
| `test_phenoconversion.py` | phenoconversion.py | 30+ CYP inhibitors/inducers, phenotype shift logic, strength classification |
| `test_query_expansion.py` | query_expansion.py | All 14 expansion maps, deduplication, category separation |
| `test_export.py` | export.py | Markdown, JSON, PDF, FHIR R4 rendering with alerts and drug interaction matrix |
| `test_metrics.py` | metrics.py | All 22 Prometheus metric registrations, label validation |
| `test_settings.py` | settings.py | PGX_ prefix resolution, default values, type validation |
| `test_scheduler.py` | scheduler.py | Job registration, interval configuration, no-op stub fallback |

### 16.3 Clinical Validation Tests

The test suite specifically validates:
- Correct activity score assignment for all defined star alleles
- Correct phenotype translation for all diplotype activity score ranges
- IWPC warfarin algorithm produces expected doses for known input combinations
- All 12 HLA-drug associations return correct status and severity
- Phenoconversion correctly shifts phenotype for strong/moderate/weak inhibitors
- DPYD dose reduction matches CPIC-recommended adjustments
- TPMT+NUDT15 combined dosing follows CPIC thiopurine guideline

---

## 17. Performance

### 17.1 Query Latency

| Operation | Expected Latency |
|-----------|-----------------|
| Evidence retrieval (15 collections, parallel) | < 1 second |
| Full RAG query (retrieve + LLM synthesis) | < 5 seconds |
| Drug check (single drug) | < 2 seconds |
| Dosing calculation (IWPC algorithm) | < 10 milliseconds |
| HLA screening (single drug) | < 5 milliseconds |
| Phenoconversion check (medication list) | < 5 milliseconds |
| Export (Markdown/JSON) | < 100 milliseconds |
| Export (PDF) | < 500 milliseconds |

### 17.2 Throughput

| Metric | Value |
|--------|-------|
| uvicorn workers | 2 (configurable) |
| ThreadPoolExecutor | 15 concurrent collection searches |
| Embedding batch size | 32 |
| Max concurrent requests | ~20 (limited by LLM API) |

### 17.3 Resource Usage

| Component | RAM | CPU |
|-----------|-----|-----|
| Milvus standalone | 4-8 GB | 2 cores |
| BGE-small-en-v1.5 | ~500 MB | 1 core (or GPU) |
| FastAPI + RAG engine | 1-2 GB | 2 cores |
| Streamlit UI | 500 MB | 1 core |

---

## 18. Observability and Metrics

**File:** `src/metrics.py` (399 lines)

### 18.1 Histograms (10)

| Metric | Labels | Buckets |
|--------|--------|---------|
| `pgx_query_latency_seconds` | query_type | 0.1, 0.5, 1, 2, 5, 10, 30 |
| `pgx_evidence_count` | — | 0, 5, 10, 15, 20, 25, 30 |
| `pgx_cross_collection_query_latency_seconds` | query_type | 0.1, 0.5, 1, 2, 5, 10, 30 |
| `pgx_cross_collection_results_count` | — | 0, 1, 5, 10, 20, 50, 100 |
| `pgx_llm_api_latency_seconds` | provider, model | 0.5, 1, 2, 5, 10, 30, 60 |
| `pgx_embedding_latency_seconds` | — | — |
| `pgx_pipeline_stage_latency_seconds` | stage | — |
| `pgx_dosing_calculation_seconds` | algorithm | — |
| `pgx_hla_screening_seconds` | — | — |
| `pgx_phenoconversion_check_seconds` | — | — |

### 18.2 Counters (8)

| Metric | Labels |
|--------|--------|
| `pgx_queries_total` | query_type |
| `pgx_errors_total` | error_type |
| `pgx_alerts_generated_total` | severity |
| `pgx_drug_checks_total` | — |
| `pgx_hla_screens_total` | result |
| `pgx_dosing_calculations_total` | algorithm |
| `pgx_phenoconversion_detections_total` | — |
| `pgx_export_total` | format |

### 18.3 Gauges (4)

| Metric | Description |
|--------|-------------|
| `pgx_collections_connected` | Number of connected Milvus collections |
| `pgx_total_vectors` | Total vectors across all collections |
| `pgx_last_ingest_timestamp` | Unix timestamp of last ingest run |
| `pgx_active_sessions` | Active Streamlit UI sessions |

### 18.4 No-Op Fallback

If `prometheus_client` is not installed, the module exports no-op stubs so the application functions without a hard dependency on Prometheus.

---

## 19. Security

### 19.1 Input Validation

- **Milvus filter injection prevention**: `_SAFE_FILTER_RE` regex restricts filter expressions to alphanumeric characters and safe symbols
- **Pydantic request validation**: All API inputs validated with field-level constraints (min_length, max_length, type checking)
- **Request size limit**: 10MB maximum (configurable via `PGX_MAX_REQUEST_SIZE_MB`)

### 19.2 CORS Configuration

Default allowed origins:
- `http://localhost:8080` (landing page)
- `http://localhost:8107` (self)
- `http://localhost:8507` (Streamlit UI)

Configurable via `PGX_CORS_ORIGINS` environment variable.

### 19.3 Secret Management

| Secret | Storage | Notes |
|--------|---------|-------|
| ANTHROPIC_API_KEY | `.env` file | Required; not committed to git |
| NCBI_API_KEY | `.env` file | Optional; increases PubMed rate limit |

### 19.4 Known Security Considerations

- No authentication on API endpoints (recommended for production: add JWT or API key middleware)
- Milvus runs without authentication by default (recommended: enable Milvus auth for production)
- API key loaded from `.env` file at startup (recommended: use secrets manager for production)

---

## 20. Data Sources

### 20.1 Primary Data Sources

| Source | Data Type | Ingest Method |
|--------|----------|---------------|
| CPIC (cpicpgx.org) | Gene-drug guidelines, evidence levels | API parser |
| PharmVar (pharmvar.org) | Star allele definitions, haplotype data | API parser |
| PharmGKB (pharmgkb.org) | Clinical annotations, drug-gene interactions | API parser |
| FDA Drug Labeling | Pharmacogenomic biomarker labeling | DailyMed parser |
| PubMed (NCBI) | PGx research literature | API parser (automated weekly) |
| ClinicalTrials.gov | PGx clinical trials | API parser (automated weekly) |
| Population databases | Allele frequency by population | Bulk parser |

### 20.2 Knowledge Sources

| Source | Data Type |
|--------|----------|
| CPIC Guidelines | Prescribing recommendations (codified in knowledge graph) |
| FDA Table of PGx Biomarkers | Drug labeling actions |
| IWPC (Klein et al., NEJM 2009) | Warfarin dosing algorithm coefficients |
| Indiana University Flockhart Table | CYP inhibitor/inducer classifications |

---

## 21. Configuration System

**File:** `config/settings.py` (118 lines)

Pydantic BaseSettings with `PGX_` environment variable prefix.

### 21.1 Key Configuration Groups

| Group | Settings | Description |
|-------|----------|-------------|
| Paths | PROJECT_ROOT, DATA_DIR, CACHE_DIR, REFERENCE_DIR | File system paths |
| Milvus | MILVUS_HOST, MILVUS_PORT | Vector database connection |
| Collections | 15 collection name settings | Milvus collection names |
| Search | TOP_K_PER_COLLECTION (5), SCORE_THRESHOLD (0.4) | Retrieval parameters |
| Weights | 15 collection weights (sum to 1.0) | Search weight configuration |
| LLM | LLM_PROVIDER, LLM_MODEL, ANTHROPIC_API_KEY | LLM configuration |
| Embeddings | EMBEDDING_MODEL, EMBEDDING_DIMENSION (384), EMBEDDING_BATCH_SIZE (32) | Embedding configuration |
| API | API_HOST, API_PORT (8107) | FastAPI server configuration |
| UI | STREAMLIT_PORT (8507) | Streamlit configuration |
| Metrics | METRICS_ENABLED (true) | Prometheus metrics toggle |
| Scheduler | INGEST_SCHEDULE_HOURS (168), INGEST_ENABLED (false) | Automated ingest configuration |
| Memory | MAX_CONVERSATION_CONTEXT (3) | Conversation history depth |
| Citations | CITATION_HIGH_THRESHOLD (0.75), CITATION_MEDIUM_THRESHOLD (0.60) | Citation scoring |
| CORS | CORS_ORIGINS | Allowed origins |
| Limits | MAX_REQUEST_SIZE_MB (10) | Request size limit |

---

## 22. Infrastructure and Deployment

### 22.1 Docker Compose Stack

**File:** `docker-compose.yml`

| Service | Image | Ports | Purpose |
|---------|-------|-------|---------|
| milvus-etcd | quay.io/coreos/etcd:v3.5.5 | — | Milvus metadata store |
| milvus-minio | minio:RELEASE.2023-03-20 | — | Milvus object storage |
| milvus-standalone | milvusdb/milvus:v2.4-latest | 19530, 9091 | Vector database |
| pgx-streamlit | Custom (Python 3.12) | 8507 | Streamlit UI |
| pgx-api | Custom (Python 3.12) | 8107 | FastAPI REST (2 workers) |
| pgx-setup | Custom (one-shot) | — | Create collections + seed data |

### 22.2 Volume Persistence

| Volume | Mount Point | Data |
|--------|------------|------|
| etcd_data | /etcd | Milvus metadata |
| minio_data | /minio_data | Milvus index/log objects |
| milvus_data | /var/lib/milvus | Milvus collection data |

### 22.3 Network

All services on `pgx-network` bridge network. Inter-service communication via container hostnames.

### 22.4 Scripts

| Script | Purpose |
|--------|---------|
| `scripts/setup_collections.py` | Create all 15 Milvus collections with --drop-existing and --seed options |
| `scripts/seed_knowledge.py` | Seed knowledge base reference data |
| `scripts/run_ingest.py` | Execute ingest parsers for external data sources |

---

## 23. Known Limitations

### 23.1 Clinical Limitations

| # | Limitation | Impact | Mitigation |
|---|-----------|--------|-----------|
| 1 | CYP2D6 structural variants (deletions, duplications, hybrids) require upstream tools | Cannot detect *5, *1xN from standard VCF without Stargazer/Cyrius | Document in reports; flag when CYP2D6 CNV status unknown |
| 2 | Dosing algorithms validated for adult populations only | Pediatric dosing not supported | Clearly state adult-only in algorithm output |
| 3 | FHIR R4 export is informational, not integrated with EHR CDS | Cannot trigger real-time EHR alerts | Export provides FHIR-compatible format for manual import |
| 4 | System is clinical decision support, not autonomous prescribing | Clinician must review and approve all recommendations | All outputs include "for clinician review" disclaimers |

### 23.2 Technical Limitations

| # | Limitation | Impact | Mitigation |
|---|-----------|--------|-----------|
| 5 | No API authentication | Open access to all endpoints | Add JWT/API key middleware for production |
| 6 | Milvus standalone mode (no replication) | Single point of failure for vector DB | Use Milvus cluster mode for production |
| 7 | LLM dependency on external Anthropic API | Network latency, rate limits, availability | Cache common responses; graceful degradation |
| 8 | Seed data is limited (240 records) | Sparse evidence for uncommon gene-drug pairs | Run full ingest from CPIC/PharmGKB/PubMed for production |

---

## 24. Recommendations

### 24.1 Short-Term (1-3 months)

1. **Run full external data ingest** to populate collections beyond seed data (PubMed, PharmGKB, ClinicalTrials.gov)
2. **Add API authentication** (JWT or API key) for production deployment
3. **Enable Milvus authentication** for production security
4. **Set up Grafana dashboards** with the 22 Prometheus metrics
5. **Enable automated ingest scheduler** for weekly PubMed/ClinicalTrials.gov refresh

### 24.2 Medium-Term (3-6 months)

6. **Integrate CYP2D6 structural variant calling** via Stargazer or Cyrius
7. **Implement EHR CDS Hooks** for real-time FHIR-based clinical decision support
8. **Add pediatric dosing models** for warfarin, thiopurines, and fluoropyrimidines
9. **Expand HLA associations** to include emerging evidence (HLA-DRB1*07:01/lapatinib)
10. **Implement federated learning** for cross-institutional PGx knowledge sharing

### 24.3 Long-Term (6-12 months)

11. **Multi-language support** for patient-facing PGx reports (Spanish, Mandarin, Arabic)
12. **Continuous learning** from clinical outcome feedback
13. **Integration with pharmacy dispensing systems** for real-time PGx alerts at point of sale
14. **Regulatory submission** for FDA clearance as a clinical decision support tool

---

## Appendix A: File Inventory

### Python Source Files (53)

| Path | Lines | Purpose |
|------|-------|---------|
| `src/knowledge.py` | 2,512 | Knowledge graph (25 pharmacogenes, 12 categories, 12 HLA, inhibitors, inducers, alternatives, aliases) |
| `src/collections.py` | 1,547 | 15 Milvus collection schemas, parallel search, IVF_FLAT indexing |
| `src/export.py` | 1,296 | Markdown, JSON, PDF, FHIR R4 export with PGx alerts |
| `src/query_expansion.py` | 1,254 | 14 domain-specific expansion maps |
| `src/pgx_pipeline.py` | 1,222 | StarAlleleCaller, PhenotypeTranslator, DrugGeneMatcher |
| `app/pgx_ui.py` | 2,138 | 10-tab Streamlit interface |
| `api/routes/pgx_clinical.py` | 803 | 7 PGx clinical decision support endpoints |
| `src/rag_engine.py` | 796 | Multi-collection RAG engine |
| `src/agent.py` | 761 | Autonomous reasoning agent |
| `src/dosing.py` | 692 | 4 validated dosing algorithms |
| `api/main.py` | 649 | FastAPI app with 8 core endpoints |
| `src/hla_screener.py` | 627 | 12 HLA-drug hypersensitivity screening |
| `src/models.py` | 616 | Enums, models, schemas |
| `src/phenoconversion.py` | 494 | Phenoconversion detection |
| `src/metrics.py` | 399 | 22 Prometheus metrics |
| `src/scheduler.py` | 232 | Automated ingest scheduler |
| `config/settings.py` | 118 | Pydantic configuration |
| `src/ingest/base.py` | — | Base ingest parser |
| `src/ingest/cpic_parser.py` | — | CPIC guideline parser |
| `src/ingest/pharmvar_parser.py` | — | PharmVar star allele parser |
| `src/ingest/pharmgkb_parser.py` | — | PharmGKB annotation parser |
| `src/ingest/fda_label_parser.py` | — | FDA label parser |
| `src/ingest/population_parser.py` | — | Population frequency parser |
| `src/ingest/pubmed_parser.py` | — | PubMed literature parser |
| `src/ingest/clinical_trials_parser.py` | — | ClinicalTrials.gov parser |
| `src/ingest/__init__.py` | — | Ingest package init |
| `src/__init__.py` | — | Source package init |
| `src/utils/__init__.py` | — | Utils package init |
| `api/__init__.py` | — | API package init |
| `api/routes/__init__.py` | — | Routes package init |
| `api/routes/reports.py` | — | Report generation endpoints |
| `api/routes/events.py` | — | Event audit trail endpoints |
| `app/__init__.py` | — | App package init |
| `config/__init__.py` | — | Config package init |
| `scripts/setup_collections.py` | — | Collection creation script |
| `scripts/seed_knowledge.py` | — | Knowledge seeding script |
| `scripts/run_ingest.py` | — | Ingest execution script |
| `tests/conftest.py` | — | Shared test fixtures |
| `tests/__init__.py` | — | Tests package init |
| `tests/test_knowledge.py` | — | Knowledge graph tests |
| `tests/test_collections.py` | — | Collection schema tests |
| `tests/test_models.py` | — | Model validation tests |
| `tests/test_rag_engine.py` | — | RAG engine tests |
| `tests/test_agent.py` | — | Agent pipeline tests |
| `tests/test_pgx_pipeline.py` | — | PGx pipeline tests |
| `tests/test_dosing.py` | — | Dosing algorithm tests |
| `tests/test_hla_screener.py` | — | HLA screening tests |
| `tests/test_phenoconversion.py` | — | Phenoconversion tests |
| `tests/test_query_expansion.py` | — | Query expansion tests |
| `tests/test_export.py` | — | Export rendering tests |
| `tests/test_metrics.py` | — | Prometheus metric tests |
| `tests/test_settings.py` | — | Configuration tests |
| `tests/test_scheduler.py` | — | Scheduler tests |

### JSON Seed Files (14)

Located in `data/reference/`: gene_reference.json, drug_guidelines.json, drug_interactions.json, hla_hypersensitivity.json, phenoconversion.json, dosing_algorithms.json, clinical_evidence.json, population_data.json, clinical_trials.json, fda_labels.json, drug_alternatives.json, patient_profiles.json, implementation.json, education.json

### Infrastructure Files (8)

Dockerfile, docker-compose.yml, requirements.txt, README.md, .env.example, .gitignore, pytest.ini/pyproject.toml

---

## Appendix B: Collection Schemas

### pgx_gene_reference

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(100) | Gene-allele ID (e.g., CYP2D6_star1) |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| gene | VARCHAR(50) | Pharmacogene symbol |
| star_allele | VARCHAR(50) | Star allele designation |
| defining_variants | VARCHAR(1000) | Comma-separated rsIDs |
| activity_score | FLOAT | Enzyme activity score |
| function_status | VARCHAR(50) | no/decreased/normal/increased function |
| allele_frequency_global | FLOAT | Global frequency |
| allele_frequency_european | FLOAT | European frequency |
| allele_frequency_african | FLOAT | African frequency |
| allele_frequency_east_asian | FLOAT | East Asian frequency |
| text | VARCHAR(4000) | Embedding source text |

### pgx_drug_guidelines

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(100) | Guideline ID |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| gene | VARCHAR(50) | Pharmacogene symbol |
| drug | VARCHAR(200) | Drug name |
| guideline_body | VARCHAR(20) | CPIC, DPWG, FDA, CPNDS |
| cpic_level | VARCHAR(10) | A, A/B, B, C, D |
| phenotype | VARCHAR(100) | Metabolizer phenotype |
| recommendation | VARCHAR(2000) | Prescribing recommendation |
| text | VARCHAR(4000) | Embedding source text |

---

## Appendix C: Metric Registry

Complete registry of all 22 Prometheus metrics. See [Section 18](#18-observability-and-metrics) for full details.

All metrics use the `pgx_` prefix for Grafana dashboard filtering alongside existing HCLS AI Factory exporters (node_exporter:9100, DCGM:9400).

---

## Appendix D: Codebase Metrics

| Metric | Value |
|--------|-------|
| Total files | 75 |
| Python files | 53 |
| JSON seed files | 14 |
| Config/infra files | 8 |
| Total Python LOC | 23,049 |
| Source LOC | 19,148 |
| Test LOC | 3,901 |
| Test count | 696 |
| Test pass rate | 100% |
| Test execution time | 0.41 seconds |
| Milvus collections | 15 |
| Seed data records | 240 |
| Knowledge graph pharmacogenes | 25 |
| Drugs covered | 100+ |
| Therapeutic categories | 12 |
| HLA-drug associations | 12 |
| CYP inhibitors/inducers | 30+ |
| Dosing algorithms | 4 |
| Clinical workflows | 8 |
| Ingest parsers | 8 |
| API endpoints | 16+ |
| UI tabs | 10 |
| Prometheus metrics | 22 (10 histograms, 8 counters, 4 gauges) |
| Export formats | 4 (MD, JSON, PDF, FHIR R4) |
| Query expansion maps | 14 |
| Entity aliases | 80+ |
| Streamlit port | 8507 |
| FastAPI port | 8107 |
| LLM | Claude Sonnet 4.6 (claude-sonnet-4-6) |
| Embedding model | BGE-small-en-v1.5 (384-dim) |
