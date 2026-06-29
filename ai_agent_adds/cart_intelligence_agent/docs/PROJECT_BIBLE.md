# CAR-T Intelligence Agent -- Project Bible

**Version:** 2.0.0
**Author:** Adam Jones (14+ years genomic research experience)
**Date:** March 2026
**License:** Apache 2.0
**Repository:** hcls-ai-factory/ai_agent_adds/cart_intelligence_agent

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Vision & Mission](#2-vision--mission)
3. [System Overview](#3-system-overview)
4. [Collections Catalog](#4-collections-catalog)
5. [Knowledge Graph](#5-knowledge-graph)
6. [Query Expansion](#6-query-expansion)
7. [Data Models](#7-data-models)
8. [RAG Engine](#8-rag-engine)
9. [Agent Architecture](#9-agent-architecture)
10. [Export System](#10-export-system)
11. [Ingest Pipelines](#11-ingest-pipelines)
12. [API Reference](#12-api-reference)
13. [UI Guide](#13-ui-guide)
14. [Metrics & Monitoring](#14-metrics--monitoring)
15. [Scheduler](#15-scheduler)
16. [Configuration](#16-configuration)
17. [Docker Deployment](#17-docker-deployment)
18. [Testing](#18-testing)
19. [Project Timeline](#19-project-timeline)
20. [File Inventory](#20-file-inventory)
21. [Dependencies](#21-dependencies)
22. [Future Roadmap](#22-future-roadmap)

---

## 1. Executive Summary

The **CAR-T Intelligence Agent** is a domain-specialized retrieval-augmented generation (RAG) system that provides cross-functional intelligence across the entire CAR-T cell therapy development lifecycle. It operates within the HCLS AI Factory Precision Intelligence Network -- one of three GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence Network, Therapeutic Discovery Engine) that together deliver patient DNA to drug candidates in under 5 hours.

The agent searches 11 Milvus vector collections containing 3,567,622 vectors -- spanning published literature, clinical trials, CAR construct designs, functional assays, manufacturing records, pharmacovigilance data, predictive biomarkers, FDA regulatory milestones, molecular sequence data, real-world evidence, and patient genomic variants -- to answer complex questions about chimeric antigen receptor T-cell therapy.

The system breaks down the data silos that typically separate target identification, CAR design, manufacturing, clinical development, and post-market surveillance. A researcher asking "Why do CD19 CAR-T therapies fail?" receives a unified answer drawing simultaneously from resistance biology literature, clinical trial outcome data, manufacturing quality records, safety databases, and genomic variant evidence.

### Platform Positioning

The HCLS AI Factory comprises three engines:

1. **Genomic Foundation Engine** -- Parabricks/DeepVariant/BWA-MEM2 (FASTQ to VCF)
2. **Precision Intelligence Network** -- 11 domain-specialized RAG agents including this CAR-T agent
3. **Therapeutic Discovery Engine** -- BioNeMo MolMIM/DiffDock/RDKit (molecular generation and docking)

### The 11 Peer Agents

| # | Agent | Domain |
|---|-------|--------|
| 1 | Precision Oncology Agent | Molecular tumor board decision support |
| 2 | **CAR-T Intelligence Agent** (this agent) | CAR-T cell therapy development lifecycle |
| 3 | Precision Biomarker Agent | Genotype-aware biomarker interpretation |
| 4 | Clinical Trial Intelligence Agent | Trial landscape monitoring and matching |
| 5 | Cardiology Intelligence Agent | Cardiovascular risk and treatment optimization |
| 6 | Neurology Intelligence Agent | Neurological disease characterization |
| 7 | Pharmacogenomics (PGx) Agent | Drug-gene interaction and dosing guidance |
| 8 | Medical Imaging Intelligence Agent | Imaging AI with NVIDIA NIM microservices |
| 9 | Single-Cell Intelligence Agent | Single-cell omics analysis |
| 10 | Autoimmune Intelligence Agent | Autoimmune disease monitoring and treatment |
| 11 | Rare Disease Intelligence Agent | Rare disease diagnosis and therapeutic matching |

### Cross-Agent Integration

The CAR-T agent calls 5 peer agents via `cross_modal/cross_agent.py` and the `/integrated-assessment` endpoint: Biomarker, Oncology, Single-Cell, Cardiology, and Clinical Trial agents.

### Pediatric CAR-T as Primary Use Case

Pediatric oncology is a primary use case for the CAR-T agent. Tisagenlecleucel (Kymriah) was the first FDA-approved CAR-T product, indicated for pediatric and young adult r/r B-ALL with 82% CR in the ELIANA trial. The agent supports pediatric-specific content including CD22 targeting for CD19-negative relapse, GD2-directed CAR-T for neuroblastoma, and pediatric CRS/ICANS management protocols with weight-based tocilizumab dosing.

The CAR-T Intelligence Agent extends the existing RAG/Chat pipeline with 10 new domain-specific collections, a 71-entry knowledge graph, 229 query expansion keywords, and an autonomous agent with plan-search-synthesize-report reasoning.

### Codebase at a Glance

| Metric | Value |
|--------|-------|
| Total Python lines | 21,259 |
| Python files | 61 |
| src/ (core engine) | 12,944 lines |
| app/ (Streamlit UI) | 1,162 lines |
| tests/ (test suite) | 4,321 lines (415 tests) |
| scripts/ (ingest/seed/setup) | 1,686 lines |
| config/ (Pydantic settings) | 113 lines |
| Milvus collections | 11 |
| Total vectors | 3,567,622 |
| Knowledge graph entries | 71 |
| Query expansion keywords | 229 |
| Query expansion terms | 1,961 |
| Seed data records | 649 |

---

## 2. Vision & Mission

### The Democratization Thesis

CAR-T cell therapy is the most successful form of adoptive cell therapy, with 6 FDA-approved products generating over $5 billion in annual revenue. Yet the intelligence required to develop a new CAR-T product is fragmented across dozens of disconnected databases, journals, regulatory filings, and institutional knowledge silos. A single researcher cannot simultaneously track target antigen biology, clinical trial outcomes, manufacturing optimization, toxicity management protocols, biomarker validation, regulatory precedent, molecular design, and real-world outcomes.

The CAR-T Intelligence Agent collapses these silos into a single conversational interface. It is designed to run on an NVIDIA DGX Spark -- a $4,699 workstation with a GB10 GPU, 128GB unified LPDDR5x memory, and 20 ARM cores -- making world-class CAR-T intelligence accessible to any lab, hospital, or biotech company.

### Mission Statement

Accelerate CAR-T cell therapy development by providing unified, evidence-grounded intelligence across the entire development lifecycle -- from target identification through post-market surveillance -- on hardware that any institution can afford.

### Design Principles

1. **Evidence-grounded**: Every answer cites specific literature, trials, assay data, or regulatory records with clickable links to source databases.
2. **Cross-functional**: Insights connect across development stages (e.g., how manufacturing choices affect clinical outcomes, how biomarkers predict safety events).
3. **Failure-mode aware**: The system proactively highlights resistance mechanisms, manufacturing failure modes, and safety signals.
4. **Quantitative**: Responses include specific numbers -- response rates, incidence percentages, binding affinities, expansion fold-changes -- not vague generalizations.
5. **Regulatory-contextualized**: Product-level intelligence always includes FDA/EMA approval status, REMS requirements, and post-marketing commitments.

---

## 3. System Overview

### High-Level Architecture

```
+------------------------------------------------------------------+
|                        STREAMLIT UI (8521)                        |
|  [Chat Tab]  [Knowledge Graph Tab]  [Image Analysis Tab]         |
|  Sidebar: Target filter, Stage filter, Date range, Collections   |
|  Demo queries | Export: MD / JSON / PDF                          |
+------------------------------------------------------------------+
         |                                          |
         v                                          v
+------------------+                    +---------------------+
| CARTIntelligence |                    | FastAPI REST (8522) |
|     Agent        |                    | 13 endpoints        |
|  plan -> search  |                    | /health /collections|
|  -> synthesize   |                    | /query /search      |
|  -> report       |                    | /find-related       |
+------------------+                    | /knowledge/stats    |
         |                              | /metrics            |
         v                              +---------------------+
+------------------------------------------------------------------+
|                      CARTRAGEngine                                |
|  embed -> parallel search -> merge & rank -> knowledge augment   |
|  -> LLM synthesis (Claude Sonnet 4.6)                            |
|                                                                  |
|  Components:                                                     |
|  +------------------+  +--------------+  +-------------------+   |
|  | QueryExpander    |  | Knowledge    |  | ExportEngine      |   |
|  | 12 maps          |  | Graph        |  | MD / JSON / PDF   |   |
|  | 229 keywords     |  | 3 dicts      |  | NVIDIA theming    |   |
|  | 1,961 terms      |  | 71 entries   |  |                   |   |
|  +------------------+  +--------------+  +-------------------+   |
+------------------------------------------------------------------+
         |
         v
+------------------------------------------------------------------+
|                 CARTCollectionManager                             |
|  11 collections | IVF_FLAT index | COSINE metric                 |
|  ThreadPoolExecutor parallel search across all collections       |
+------------------------------------------------------------------+
         |
         v
+------------------------------------------------------------------+
|                    Milvus 2.4 (19530)                            |
|  BGE-small-en-v1.5 embeddings (384-dim)                          |
|                                                                  |
|  cart_literature    (5,047)    cart_safety      (71)             |
|  cart_trials        (973)     cart_biomarkers  (60)             |
|  cart_constructs    (41)      cart_regulatory  (40)             |
|  cart_assays        (75)      cart_sequences   (40)             |
|  cart_manufacturing (56)      cart_realworld   (54)             |
|                                                                  |
|  genomic_evidence   (3,561,170) [read-only, shared]              |
+------------------------------------------------------------------+
         |
         v
+------------------------------------------------------------------+
|                   Claude Sonnet 4.6 (Anthropic API)              |
|  System prompt: 12 expertise domains                             |
|  Streaming response | Citation-linked output                     |
+------------------------------------------------------------------+
```

### Data Flow

```
User Question
     |
     v
[1] Search Plan (agent identifies targets, stages, strategy)
     |
     v
[2] Query Embedding (BGE-small-en-v1.5, "Represent this sentence...")
     |
     v
[3] Parallel Search (ThreadPoolExecutor across 11 collections)
     |
     v
[4] Query Expansion (12 maps, re-embed top 5 expansion terms)
     |
     v
[5] Merge & Rank (deduplicate, weighted scoring, cap at 30 hits)
     |
     v
[6] Knowledge Augmentation (targets, toxicities, manufacturing,
     |  biomarkers, regulatory, immunogenicity)
     v
[7] Prompt Assembly (evidence sections + knowledge context + question)
     |
     v
[8] LLM Synthesis (Claude Sonnet 4.6, streaming, 2048 tokens)
     |
     v
[9] Response with clickable citations + evidence panel + export buttons
```

---

## 4. Collections Catalog

The system manages 11 Milvus collections: 10 owned by the CAR-T Intelligence Agent and 1 shared read-only collection from the upstream rag-chat-pipeline.

### 4.1 Collection Summary

| # | Collection | Vectors | Description | Data Sources |
|---|-----------|---------|-------------|-------------|
| 1 | `cart_literature` | 5,047 | Published research papers and patents | PubMed, PMC, patents |
| 2 | `cart_trials` | 973 | Clinical trial records | ClinicalTrials.gov |
| 3 | `cart_constructs` | 41 | CAR construct designs and engineering approaches | Seed data |
| 4 | `cart_assays` | 75 | In vitro/in vivo functional assay results | Seed data |
| 5 | `cart_manufacturing` | 56 | CMC/manufacturing process records | Seed data |
| 6 | `cart_safety` | 71 | Adverse events, CRS/ICANS incidence | Seed data, FAERS |
| 7 | `cart_biomarkers` | 60 | Predictive/prognostic/PD biomarkers | Seed data |
| 8 | `cart_regulatory` | 40 | FDA approvals, EMA, PMDA, NMPA designations | Seed data |
| 9 | `cart_sequences` | 40 | scFv/CAR molecular & structural data | Seed data |
| 10 | `cart_realworld` | 54 | CIBMTR, post-market outcomes | Seed data |
| 11 | `genomic_evidence` | 3,561,170 | Patient variant data (ClinVar + AlphaMissense) | Read-only, rag-chat-pipeline |
| | **TOTAL** | **3,567,622** | | |

### 4.2 Shared Configuration

All collections share:
- **Embedding model**: BGE-small-en-v1.5 (384-dimensional)
- **Index type**: IVF_FLAT
- **Index params**: `nlist=1024`
- **Search params**: `nprobe=16`
- **Metric**: COSINE similarity
- **Primary key**: VARCHAR `id` field

### 4.3 Collection Schemas (Detail)

#### cart_literature

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | PMID or patent number |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| title | VARCHAR | 500 | Paper or patent title |
| text_chunk | VARCHAR | 3000 | Text chunk for embedding |
| source_type | VARCHAR | 20 | pubmed, pmc, patent, preprint, manual |
| year | INT64 | -- | Publication year |
| cart_stage | VARCHAR | 30 | Development stage |
| target_antigen | VARCHAR | 100 | e.g., CD19, BCMA |
| disease | VARCHAR | 200 | Disease or indication |
| keywords | VARCHAR | 1000 | Comma-separated keywords/MeSH |
| journal | VARCHAR | 200 | Journal name |

#### cart_trials

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 20 | NCT number (e.g., NCT03958656) |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| title | VARCHAR | 500 | Official trial title |
| text_summary | VARCHAR | 3000 | Brief summary for embedding |
| phase | VARCHAR | 30 | Trial phase |
| status | VARCHAR | 30 | Recruitment status |
| sponsor | VARCHAR | 200 | Lead sponsor |
| target_antigen | VARCHAR | 100 | Target antigen |
| car_generation | VARCHAR | 20 | 1st, 2nd, 3rd, 4th, armored, universal |
| costimulatory | VARCHAR | 50 | CD28, 4-1BB, dual |
| disease | VARCHAR | 200 | Disease or indication |
| enrollment | INT64 | -- | Target enrollment count |
| start_year | INT64 | -- | Study start year |
| outcome_summary | VARCHAR | 2000 | Outcome summary |

#### cart_constructs

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | Construct identifier |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| name | VARCHAR | 200 | Product or construct name |
| text_summary | VARCHAR | 2000 | Description |
| target_antigen | VARCHAR | 100 | Target antigen |
| scfv_origin | VARCHAR | 200 | scFv antibody clone/origin |
| costimulatory_domain | VARCHAR | 100 | Costimulatory domain(s) |
| signaling_domain | VARCHAR | 100 | Default: CD3-zeta |
| generation | VARCHAR | 20 | CAR generation |
| hinge_tm | VARCHAR | 200 | Hinge + transmembrane domain |
| vector_type | VARCHAR | 50 | Lentiviral, retroviral, etc. |
| fda_status | VARCHAR | 20 | FDA approval status |
| known_toxicities | VARCHAR | 500 | CRS, ICANS, etc. |

#### cart_assays

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | Assay record identifier |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| text_summary | VARCHAR | 2000 | Assay description |
| assay_type | VARCHAR | 30 | cytotoxicity, cytokine, flow, etc. |
| construct_id | VARCHAR | 100 | FK to cart_constructs |
| target_antigen | VARCHAR | 100 | Target antigen tested |
| cell_line | VARCHAR | 100 | e.g., Nalm-6, Raji, K562 |
| effector_ratio | VARCHAR | 20 | E:T ratio |
| key_metric | VARCHAR | 50 | Primary metric name |
| metric_value | FLOAT | -- | Numeric value |
| outcome | VARCHAR | 20 | success, partial, failure |
| notes | VARCHAR | 1000 | Additional notes |

#### cart_manufacturing

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | Manufacturing record ID |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| text_summary | VARCHAR | 2000 | Process description |
| process_step | VARCHAR | 30 | transduction, expansion, etc. |
| vector_type | VARCHAR | 50 | Vector type |
| parameter | VARCHAR | 100 | Process parameter name |
| parameter_value | VARCHAR | 50 | Parameter value |
| target_spec | VARCHAR | 100 | Acceptance criteria |
| met_spec | VARCHAR | 10 | yes, no, borderline |
| batch_id | VARCHAR | 50 | Batch identifier |
| notes | VARCHAR | 1000 | Additional notes |

#### cart_safety

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | Safety record ID |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| text_summary | VARCHAR | 3000 | Event description |
| product | VARCHAR | 200 | Product name |
| event_type | VARCHAR | 30 | CRS, ICANS, cytopenia, etc. |
| severity_grade | VARCHAR | 100 | Grade 1-5 or mild/moderate/severe |
| onset_timing | VARCHAR | 100 | e.g., median day 5 |
| incidence_rate | VARCHAR | 200 | e.g., 42% any grade |
| management_protocol | VARCHAR | 500 | Treatment protocol |
| outcome | VARCHAR | 100 | Outcome description |
| reporting_source | VARCHAR | 50 | FAERS, trial, registry, label |
| year | INT64 | -- | Reporting year |

#### cart_biomarkers

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | Biomarker record ID |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| text_summary | VARCHAR | 3000 | Biomarker description |
| biomarker_name | VARCHAR | 100 | Biomarker name |
| biomarker_type | VARCHAR | 30 | predictive, prognostic, etc. |
| assay_method | VARCHAR | 100 | ELISA, flow cytometry, etc. |
| clinical_cutoff | VARCHAR | 100 | e.g., >500 mg/L |
| predictive_value | VARCHAR | 200 | Clinical correlation |
| associated_outcome | VARCHAR | 200 | Outcome it predicts |
| target_antigen | VARCHAR | 100 | Target antigen |
| disease | VARCHAR | 200 | Disease context |
| evidence_level | VARCHAR | 20 | validated, emerging, exploratory |

#### cart_regulatory

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | Regulatory record ID |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| text_summary | VARCHAR | 3000 | Milestone description |
| product | VARCHAR | 200 | Product name |
| regulatory_event | VARCHAR | 50 | BLA, breakthrough therapy, etc. |
| date | VARCHAR | 20 | YYYY-MM-DD or YYYY-MM |
| agency | VARCHAR | 20 | FDA, EMA, PMDA |
| indication | VARCHAR | 200 | Approved indication |
| decision | VARCHAR | 100 | approved, rejected, pending |
| conditions | VARCHAR | 500 | Approval conditions |
| pivotal_trial | VARCHAR | 100 | Pivotal trial ID |

#### cart_sequences

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | Sequence record ID |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| text_summary | VARCHAR | 3000 | Molecular description |
| construct_name | VARCHAR | 200 | Construct name |
| target_antigen | VARCHAR | 100 | Target antigen |
| scfv_clone | VARCHAR | 100 | e.g., FMC63, SJ25C1 |
| binding_affinity_kd | VARCHAR | 50 | e.g., 0.3 nM |
| heavy_chain_vregion | VARCHAR | 500 | VH framework/CDR info |
| light_chain_vregion | VARCHAR | 500 | VL framework/CDR info |
| framework | VARCHAR | 100 | IgG1, IgG4, etc. |
| species_origin | VARCHAR | 30 | murine, humanized, fully_human |
| immunogenicity_risk | VARCHAR | 20 | low, moderate, high |
| structural_notes | VARCHAR | 1000 | Structural notes |

#### cart_realworld

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 100 | RWE record ID |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| text_summary | VARCHAR | 3000 | Study description |
| study_type | VARCHAR | 30 | retrospective, registry, etc. |
| data_source | VARCHAR | 100 | CIBMTR, institutional, SEER |
| product | VARCHAR | 200 | Product name |
| indication | VARCHAR | 200 | Disease indication |
| population_size | INT64 | -- | Study population |
| median_followup_months | FLOAT | -- | Median follow-up |
| primary_endpoint | VARCHAR | 100 | Primary endpoint |
| outcome_value | VARCHAR | 100 | Outcome value |
| setting | VARCHAR | 50 | academic, community, both |
| special_population | VARCHAR | 200 | elderly, bridging, CNS, etc. |

#### genomic_evidence (Read-Only)

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 200 | Variant identifier |
| embedding | FLOAT_VECTOR | 384 | BGE-small-en-v1.5 |
| chrom | VARCHAR | 10 | Chromosome |
| pos | INT64 | -- | Position |
| ref | VARCHAR | 500 | Reference allele |
| alt | VARCHAR | 500 | Alternate allele |
| qual | FLOAT | -- | Quality score |
| gene | VARCHAR | 50 | Gene symbol |
| consequence | VARCHAR | 100 | Variant consequence |
| impact | VARCHAR | 20 | HIGH, MODERATE, etc. |
| genotype | VARCHAR | 10 | e.g., 0/1, 1/1 |
| text_summary | VARCHAR | 2000 | Summary text |
| clinical_significance | VARCHAR | 200 | ClinVar significance |
| rsid | VARCHAR | 20 | dbSNP rsID |
| disease_associations | VARCHAR | 500 | Disease associations |
| am_pathogenicity | FLOAT | -- | AlphaMissense score |
| am_class | VARCHAR | 30 | AlphaMissense class |

---

## 5. Knowledge Graph

The knowledge graph provides structured, curated domain knowledge that augments vector search results before LLM synthesis. It contains 3 dictionaries with 71 total entries and 54 entity aliases for cross-entity resolution.

### 5.1 Dictionary Summary

| Dictionary | Entries | Description |
|-----------|---------|-------------|
| CART_TARGETS | 34 | Target antigen profiles with expression, diseases, products, resistance |
| CART_TOXICITIES | 17 | Toxicity profiles with grading, management, biomarkers, risk factors |
| CART_MANUFACTURING | 20 | Manufacturing processes with parameters, failure modes, specs |
| CART_BIOMARKERS | 23 | Predictive/prognostic biomarkers with cutoffs and evidence levels |
| CART_REGULATORY | 6 | FDA-approved products with approval history and REMS |
| CART_IMMUNOGENICITY | 6 | HLA, ADA, humanization, and immunogenicity testing |
| **TOTAL** | **71** | |

### 5.2 CART_TARGETS (34 Entries)

Each entry includes: protein name, UniProt ID, expression pattern, diseases, approved products, key trials, resistance mechanisms, toxicity profile, and normal tissue expression.

**Targets with approved products (2):**
- **CD19**: Kymriah, Yescarta, Tecartus, Breyanzi (B-ALL, DLBCL, FL, MCL, CLL)
- **BCMA**: Abecma, Carvykti (Multiple Myeloma)

**Hematologic targets (13):** CD19, BCMA, CD22, CD20, CD30, CD33, CD38, CD123, GPRC5D, NKG2D_ligands, CD7, CD5, ROR1

**Solid tumor targets (12):** GD2, HER2, GPC3, EGFR, EGFRvIII, Mesothelin, Claudin18.2, MUC1, PSMA, IL13Ra2, DLL3, B7-H3

### 5.3 CART_TOXICITIES (17 Entries)

| Toxicity | Full Name | Incidence | Key Management |
|----------|-----------|-----------|----------------|
| CRS | Cytokine Release Syndrome | 50-95% any grade; 10-25% grade 3+ | Tocilizumab, corticosteroids |
| ICANS | Immune Effector Cell-Associated Neurotoxicity | 20-65% any grade; 10-30% grade 3+ | Corticosteroids, ICU monitoring |
| B_CELL_APLASIA | B-Cell Aplasia / Hypogammaglobulinemia | Near 100% with CD19 CAR-T | IVIG replacement |
| HLH_MAS | Hemophagocytic Lymphohistiocytosis / MAS | 1-5% | Etoposide, anakinra, ruxolitinib |
| CYTOPENIAS | Prolonged Cytopenias | 30-50% prolonged (>30d) | G-CSF, transfusion support |
| TLS | Tumor Lysis Syndrome | <5% | Rasburicase, hydration |
| GVHD | Graft-versus-Host Disease | Rare autologous; risk allogeneic | TCR KO prevention, steroids |
| ON_TARGET_OFF_TUMOR | On-Target, Off-Tumor Toxicity | Target-dependent | Affinity tuning, logic gates, safety switches |

### 5.4 CART_MANUFACTURING (20 Entries)

lentiviral_transduction, retroviral_transduction, t_cell_activation, ex_vivo_expansion, leukapheresis, cryopreservation, release_testing, point_of_care_manufacturing, lymphodepletion, vein_to_vein_time

### 5.5 CART_BIOMARKERS (23 Entries)

| Biomarker | Type | Clinical Cutoff | Evidence Level |
|-----------|------|----------------|----------------|
| Ferritin | Predictive | >500 mg/L pre-infusion | Validated |
| CRP | Predictive | >200 mg/L within 72h | Validated |
| IL-6 | Pharmacodynamic | >1000 pg/mL | Validated |
| sIL-2R | Pharmacodynamic | >5000 pg/mL | Validated |
| CAR-T Expansion (Cmax) | Pharmacodynamic | >50,000 copies/ug DNA | Validated |
| Tcm% | Predictive | >40% in apheresis | Emerging |
| CD4:CD8 Ratio | Predictive | Optimal 1:1 | Emerging |
| LDH | Prognostic | >ULN | Validated |
| PD-1 | Resistance | >30% PD-1+ on CAR-T | Emerging |
| LAG-3 | Resistance | >20% LAG-3+ | Emerging |
| TIM-3 | Resistance | >25% TIM-3+ | Emerging |
| MRD (flow) | Monitoring | <10^-4 (0.01%) | Validated |
| ctDNA | Monitoring | 2-log reduction by day 28 | Emerging |
| sBCMA | Resistance | >40 ng/mL baseline | Emerging |
| IFN-gamma | Pharmacodynamic | >500 pg/mL at day 7 | Validated |

### 5.6 CART_REGULATORY (6 FDA-Approved Products)

| Product | Generic Name | Manufacturer | Initial Approval | Target |
|---------|-------------|-------------|-----------------|--------|
| Kymriah | tisagenlecleucel | Novartis | 2017-08-30 | CD19 |
| Yescarta | axicabtagene ciloleucel | Kite/Gilead | 2017-10-18 | CD19 |
| Tecartus | brexucabtagene autoleucel | Kite/Gilead | 2020-07-24 | CD19 |
| Breyanzi | lisocabtagene maraleucel | BMS/Juno | 2021-02-05 | CD19 |
| Abecma | idecabtagene vicleucel | BMS/bluebird | 2021-03-26 | BCMA |
| Carvykti | ciltacabtagene autoleucel | Janssen/Legend | 2022-02-28 | BCMA |

### 5.7 CART_IMMUNOGENICITY (6 Entries)

murine_scfv_immunogenicity, humanization_strategies, ada_clinical_impact, hla_restricted_epitopes, immunogenicity_testing, allogeneic_hla_considerations

### 5.8 Entity Aliases (ENTITY_ALIASES)

54 aliases map common names to canonical entities for cross-entity resolution in comparative queries. Categories include:

- **Product aliases** (12): Kymriah/tisagenlecleucel, Yescarta/axicabtagene, etc.
- **Costimulatory domains** (4): 4-1BB/CD137, CD28, OX40, ICOS
- **Vector types** (2): LENTIVIRAL, RETROVIRAL
- **Biomarker aliases** (13): FERRITIN, CRP, IL-6, PD-1, LAG-3, TIM-3, MRD, etc.
- **Immunogenicity aliases** (8): ADA, HAMA, HUMANIZATION, HLA, MHC, etc.

### 5.9 Context Retrieval Functions

| Function | Input | Output |
|----------|-------|--------|
| `get_target_context(antigen)` | Target name (e.g., "CD19") | Formatted markdown with biology, diseases, products, resistance |
| `get_toxicity_context(toxicity)` | Toxicity ID (e.g., "CRS") | Mechanism, grading, management, biomarkers, risk factors |
| `get_manufacturing_context(process)` | Process ID | Description, parameters, failure modes, specs |
| `get_biomarker_context(biomarker)` | Biomarker key | Type, assay, cutoff, predictive value, evidence level |
| `get_regulatory_context(product)` | Product name | Approval history, designations, REMS, EMA |
| `get_immunogenicity_context(topic)` | Topic key | Description, methods, clinical impact |
| `get_all_context_for_query(query)` | Free-text query | Combined context from all 6 domains |
| `resolve_comparison_entity(text)` | Entity name | Dict with type, canonical name, target |
| `get_comparison_context(a, b)` | Two entity dicts | Side-by-side knowledge for comparative analysis |
| `get_knowledge_stats()` | None | Dict with entry counts per dictionary |

---

## 6. Query Expansion

The query expansion system maps user keywords to related biomedical terms, broadening search recall across collections. It contains 12 expansion maps with a total of 229 keywords expanding to 1,961 terms.

### 6.1 Expansion Map Summary

| # | Map Name | Keywords | Total Terms | Domain |
|---|---------|----------|-------------|--------|
| 1 | Target Antigen | 32 | 244 | Antigen names, diseases, products, aliases |
| 2 | Disease | 16 | 143 | Disease names, related antigens, therapies |
| 3 | Toxicity | 18 | 168 | Safety events, management, biomarkers |
| 4 | Manufacturing | 21 | 215 | CMC terms, platforms, release testing |
| 5 | Mechanism | 19 | 224 | Resistance, exhaustion, signaling, TME |
| 6 | Construct | 20 | 206 | scFv, CAR generations, safety switches |
| 7 | Safety | 8 | 69 | Pharmacovigilance, AEs, REMS |
| 8 | Biomarker | 18 | 117 | Predictive markers, exhaustion, MRD |
| 9 | Regulatory | 8 | 46 | FDA, BLA, RMAT, EMA pathways |
| 10 | Sequence | 8 | 54 | scFv, CDR, binding affinity, nanobody |
| 11 | RealWorld | 10 | 65 | RWE, registries, disparities, follow-up |
| 12 | Immunogenicity | 12 | 98 | HLA, ADA, humanization, ELISpot |
| | **TOTAL** | **229** (unique) | **1,961** (unique) | |

### 6.2 Expansion Logic

The `expand_query()` function:

1. Converts the user query to lowercase
2. Iterates through all 12 expansion maps in order
3. For each keyword that matches as a substring of the query, adds all associated terms to a set
4. Returns a deduplicated, sorted list of expansion terms

The `expand_query_by_category()` variant groups results by category, enabling collection-specific weighting.

### 6.3 How Expansion Feeds Search

In the RAG engine, expansion terms are used in two ways:

1. **Antigen field filtering**: If an expansion term matches a known target antigen (from `_KNOWN_ANTIGENS`), it is used as a Milvus `target_antigen == "..."` field filter on collections that have that field.
2. **Semantic re-embedding**: Non-antigen expansion terms are re-embedded and used for additional vector similarity searches across all collections, with a 0.7x score discount to avoid drowning out primary results.

---

## 7. Data Models

The data model layer uses Pydantic v2 with strict validation. All models are defined in `src/models.py`.

### 7.1 Enums (14)

| Enum | Values | Used By |
|------|--------|---------|
| `CARTStage` | target_id, car_design, vector_eng, testing, clinical | CARTLiterature, AgentQuery |
| `SourceType` | pubmed, pmc, patent, preprint, manual | CARTLiterature |
| `TrialPhase` | Early Phase 1, Phase 1, Phase 1/Phase 2, Phase 2, Phase 2/Phase 3, Phase 3, Phase 4, N/A | ClinicalTrial |
| `TrialStatus` | Recruiting, Active not recruiting, Completed, Terminated, Withdrawn, Suspended, Not yet recruiting, Unknown | ClinicalTrial |
| `CARGeneration` | 1st, 2nd, 3rd, 4th, armored, universal | ClinicalTrial, CARConstruct |
| `AssayType` | cytotoxicity, cytokine, flow, proliferation, in_vivo, persistence, exhaustion, migration, trafficking, serial_killing | AssayResult |
| `ProcessStep` | transduction, expansion, harvest, formulation, release, cryo, non_viral, mrna_electroporation, crispr_knock_in, ipsc_derived, automated | ManufacturingRecord |
| `FDAStatus` | approved, bla_filed, phase3, phase2, phase1, preclinical, discontinued | CARConstruct |
| `SafetyEventType` | CRS, ICANS, cytopenia, infection, secondary_malignancy, organ_toxicity, neurologic, cardiac, coagulopathy, renal | SafetyRecord |
| `BiomarkerType` | predictive, prognostic, pharmacodynamic, monitoring, resistance | BiomarkerRecord |
| `EvidenceLevel` | validated, emerging, exploratory | BiomarkerRecord |
| `RegulatoryEvent` | BLA, breakthrough_therapy, RMAT, accelerated_approval, full_approval, label_update, REMS, post_marketing_requirement, complete_response | RegulatoryRecord |
| `RWEStudyType` | retrospective, registry, claims, ehr_analysis, meta_analysis | RealWorldRecord |

### 7.2 Collection Models (10)

Each collection model is a Pydantic `BaseModel` with:
- Typed, validated fields mapping to Milvus schema columns
- A `to_embedding_text() -> str` method that assembles the text passed to BGE-small-en-v1.5

| Model | Collection | Key Fields |
|-------|-----------|------------|
| `CARTLiterature` | cart_literature | id, title, text_chunk, source_type, year, cart_stage, target_antigen |
| `ClinicalTrial` | cart_trials | id (NCT), title, text_summary, phase, status, sponsor, enrollment |
| `CARConstruct` | cart_constructs | id, name, target_antigen, scfv_origin, costimulatory_domain, generation |
| `AssayResult` | cart_assays | id, text_summary, assay_type, cell_line, key_metric, metric_value, outcome |
| `ManufacturingRecord` | cart_manufacturing | id, text_summary, process_step, parameter, target_spec, met_spec |
| `SafetyRecord` | cart_safety | id, text_summary, product, event_type, severity_grade, incidence_rate |
| `BiomarkerRecord` | cart_biomarkers | id, text_summary, biomarker_name, biomarker_type, clinical_cutoff |
| `RegulatoryRecord` | cart_regulatory | id, text_summary, product, regulatory_event, date, agency, decision |
| `SequenceRecord` | cart_sequences | id, text_summary, scfv_clone, binding_affinity_kd, species_origin |
| `RealWorldRecord` | cart_realworld | id, text_summary, study_type, data_source, population_size, primary_endpoint |

### 7.3 Search & Agent Models (4)

| Model | Purpose | Key Fields |
|-------|---------|------------|
| `AgentQuery` | Input to agent/engine | question, target_antigen, cart_stage, include_genomic |
| `SearchHit` | Single search result | collection, id, score, text, metadata |
| `CrossCollectionResult` | Merged multi-collection results | query, hits[], knowledge_context, total_collections_searched, search_time_ms |
| `AgentResponse` | Full agent output | question, answer, evidence (CrossCollectionResult), knowledge_used[], timestamp |

Additionally, `ComparativeResult` extends the model set for head-to-head analysis, containing `entity_a`, `entity_b`, `evidence_a`, `evidence_b`, and `comparison_context`.

---

## 8. RAG Engine

The `CARTRAGEngine` class (defined in `src/rag_engine.py`) is the central orchestrator for all retrieval and synthesis operations.

### 8.1 Constructor

```python
CARTRAGEngine(
    collection_manager: CARTCollectionManager,
    embedder,          # SentenceTransformer wrapper, embed_text() -> List[float]
    llm_client,        # Anthropic client wrapper, generate() / generate_stream()
    knowledge,         # src.knowledge module
    query_expander,    # src.query_expansion module
)
```

### 8.2 Collection Weights

Each collection has a configurable weight that boosts its score contribution. Weights are defined in `CARTSettings` and loaded at engine initialization:

| Collection | Weight | Label |
|-----------|--------|-------|
| cart_literature | 0.20 | Literature |
| cart_trials | 0.16 | Trial |
| cart_constructs | 0.10 | Construct |
| cart_assays | 0.09 | Assay |
| cart_manufacturing | 0.07 | Manufacturing |
| cart_safety | 0.08 | Safety |
| cart_biomarkers | 0.08 | Biomarker |
| cart_regulatory | 0.06 | Regulatory |
| cart_sequences | 0.06 | Sequence |
| cart_realworld | 0.07 | RealWorld |
| genomic_evidence | 0.04 | Genomic |

The weighted score formula is: `final_score = raw_cosine_similarity * (1 + weight)`

### 8.3 Retrieval Pipeline

The `retrieve()` method executes the following steps:

1. **Query preparation**: Optionally prepend conversation context for follow-up queries.
2. **Embedding**: Embed query with BGE instruction prefix: `"Represent this sentence for searching relevant passages: "`
3. **Collection selection**: Use caller-specified filter or search all 11 collections.
4. **Filter construction**: Build per-collection Milvus filter expressions (target_antigen equality, year range bounds).
5. **Parallel search**: `CARTCollectionManager.search_all()` uses `ThreadPoolExecutor` to search all collections concurrently.
6. **Query expansion**: Call `expand_query()`, then for the top 5 expansion terms, either (a) use as field filter for known antigens or (b) re-embed and search all collections for semantic terms.
7. **Merge & rank**: Deduplicate by ID, apply weighted scoring, apply citation relevance labels (high >= 0.75, medium >= 0.60, low < 0.60), sort descending, cap at 30.
8. **Knowledge augmentation**: Scan query for target antigens, toxicities, manufacturing terms, biomarkers, products, and immunogenicity keywords; retrieve structured context from the knowledge graph.
9. **Return**: `CrossCollectionResult` with hits, knowledge context, collection count, and timing.

### 8.4 Citation Formatting

Citations are formatted as clickable Markdown links:
- **Literature**: `[Literature:PMID 12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/)`
- **Trials**: `[Trial:NCT12345678](https://clinicaltrials.gov/study/NCT12345678)`
- **Other collections**: `[Collection:record-id]`

### 8.5 Comparative Analysis

When a query contains "compare", "vs", or "versus", the engine:

1. Parses two entities from the question using regex.
2. Resolves each entity against the knowledge graph via `resolve_comparison_entity()`.
3. Runs two independent retrieval passes (one per entity).
4. Builds side-by-side knowledge context via `get_comparison_context()`.
5. Constructs a specialized comparative prompt requiring: comparison table, advantages, limitations, and clinical context.

### 8.6 System Prompt

The system prompt defines 12 expertise domains:

1. Target Identification
2. CAR Design
3. Vector Engineering
4. In Vitro & In Vivo Testing
5. Clinical Development
6. Manufacturing
7. Safety & Pharmacovigilance
8. Biomarkers
9. Regulatory Intelligence
10. Molecular Design
11. Real-World Evidence
12. Genomic Evidence

The prompt instructs the LLM to cite evidence with clickable links, think cross-functionally, highlight failure modes, be specific with quantitative data, include regulatory context, and acknowledge uncertainty.

### 8.7 Streaming

The `query_stream()` method yields three message types:
1. `{"type": "evidence", "content": CrossCollectionResult}` -- the evidence payload
2. `{"type": "token", "content": "..."}` -- individual LLM output tokens
3. `{"type": "done", "content": "full answer"}` -- the complete response

---

## 9. Agent Architecture

The `CARTIntelligenceAgent` class (defined in `src/agent.py`) wraps the RAG engine with autonomous planning and reasoning.

### 9.1 SearchPlan

```python
@dataclass
class SearchPlan:
    question: str
    identified_topics: List[str]
    target_antigens: List[str]       # Extracted from question
    relevant_stages: List[CARTStage] # Mapped from keywords
    search_strategy: str             # "broad", "targeted", or "comparative"
    sub_questions: List[str]         # Decomposed for complex queries
```

### 9.2 Agent Pipeline

```
[1] PLAN
    - Extract target antigens from question (34 known antigens)
    - Map keywords to CAR-T development stages (5 stages, ~4 keywords each)
    - Determine strategy: comparative (vs/compare), targeted (single antigen + stage), broad
    - Decompose complex questions ("why fail?") into sub-questions

[2] SEARCH
    - Build AgentQuery with extracted target_antigen
    - Call rag_engine.retrieve() with full parameters

[3] EVALUATE
    - "sufficient": >= 3 collections with hits AND >= 10 total hits
    - "partial": >= 2 collections AND >= 5 hits
    - "insufficient": anything less

[4] EXPAND (if insufficient)
    - Execute up to 2 sub-questions
    - Merge additional hits into evidence

[5] SYNTHESIZE
    - Call rag_engine.query() for LLM generation

[6] REPORT
    - Build AgentResponse with answer, evidence, and knowledge_used
    - Optionally generate a structured markdown report via generate_report()
```

### 9.3 Report Generation

`generate_report()` produces a structured Markdown document with:
- Query metadata (timestamp, collections searched, evidence count, search time)
- Full analysis text
- Evidence sources grouped by collection (top 5 per collection)
- Knowledge graph entries used

---

## 10. Export System

The export system (defined in `src/export.py`) provides three output formats with consistent structure.

### 10.1 Markdown Export

`export_markdown()` generates a complete report with:
- Title, query, timestamp, filters
- Full response text
- Evidence tables (collection-specific columns)
- Knowledge graph context
- Search metrics table
- Footer with version

### 10.2 JSON Export

`export_json()` generates structured JSON via Pydantic `.model_dump()`:
- Report type, version, timestamp
- Query and response text
- Full evidence payload (or comparative payload with entity_a/entity_b)
- Search metrics
- Applied filters

### 10.3 PDF Export

`export_pdf()` uses reportlab Platypus with NVIDIA theming:

**Color Palette:**
- NVIDIA Green: `#76B900` (headings, table headers, accent elements)
- Dark Background: `#1a1a1a`
- Table Alternating: `#f0f0f0`
- Light Gray: `#666666` (metadata, footer)

**Features:**
- Collection-specific evidence tables with alternating row colors
- Clickable PubMed and ClinicalTrials.gov links in PDF
- Markdown-to-PDF conversion (headings, bold, bullet lists, blockquotes, tables)
- Inline comparison table rendering for comparative queries
- Search metrics summary table

**Filenames:** `cart_query_YYYYMMDD_HHMMSS.{md,json,pdf}`

---

## 11. Ingest Pipelines

### 11.1 BaseIngestPipeline Pattern

All ingest pipelines inherit from `BaseIngestPipeline` (defined in `src/ingest/base.py`), which enforces a three-step workflow:

```
fetch(**kwargs)      -> Raw data (XML, JSON, CSV)
parse(raw_data)      -> List[PydanticModel]
embed_and_store()    -> Milvus insertion
```

The `embed_and_store()` method:
1. Calls each record's `to_embedding_text()` method
2. Batch-encodes texts with the embedder (batch_size=32)
3. Converts Enum values to strings
4. Truncates strings to safe UTF-8 byte lengths for Milvus VARCHAR limits
5. Inserts batch into the target collection via `CARTCollectionManager.insert_batch()`

The `run()` method orchestrates the full pipeline: `fetch -> parse -> embed_and_store`.

### 11.2 Implemented Parsers (14)

| Parser | Module | Target Collection | Data Source |
|--------|--------|------------------|-------------|
| `PubMedIngestPipeline` | literature_parser.py | cart_literature | NCBI E-utilities (PubMed API) |
| `ClinicalTrialsIngestPipeline` | clinical_trials_parser.py | cart_trials | ClinicalTrials.gov v2 API |
| `ConstructIngestPipeline` | construct_parser.py | cart_constructs | Manual curation |
| `AssayIngestPipeline` | assay_parser.py | cart_assays | JSON seed data |
| `ManufacturingIngestPipeline` | manufacturing_parser.py | cart_manufacturing | JSON seed data |
| `SafetyIngestPipeline` | safety_parser.py | cart_safety | JSON seed data |
| `BiomarkerIngestPipeline` | biomarker_parser.py | cart_biomarkers | JSON seed data |
| `RegulatoryIngestPipeline` | regulatory_parser.py | cart_regulatory | JSON seed data |
| `SequenceIngestPipeline` | sequence_parser.py | cart_sequences | JSON seed data |
| `RealWorldIngestPipeline` | realworld_parser.py | cart_realworld | JSON seed data |
| `FAERSIngestPipeline` | faers_parser.py | cart_safety | FDA FAERS database |
| `DailyMedIngestPipeline` | dailymed_parser.py | cart_regulatory | DailyMed/SPL labels |
| `UniProtIngestPipeline` | uniprot_parser.py | cart_sequences | UniProt REST API |
| `CIBMTRIngestPipeline` | cibmtr_parser.py | cart_realworld | CIBMTR registry data |

### 11.3 Seed Scripts (14)

| Script | Records | Target Collection |
|--------|---------|------------------|
| `scripts/setup_collections.py` | -- | Creates all 11 collections, optionally seeds constructs |
| `scripts/seed_knowledge.py` | -- | Validates knowledge graph integrity |
| `scripts/seed_assays.py` | 75 | cart_assays |
| `scripts/seed_manufacturing.py` | 56 | cart_manufacturing |
| `scripts/seed_safety.py` | 71 | cart_safety |
| `scripts/seed_biomarkers.py` | 60 | cart_biomarkers |
| `scripts/seed_regulatory.py` | 40 | cart_regulatory |
| `scripts/seed_sequences.py` | 40 | cart_sequences |
| `scripts/seed_realworld.py` | 54 | cart_realworld |
| `scripts/seed_patents.py` | 45 | cart_literature (patents) |
| `scripts/seed_immunogenicity.py` | 38 (20+18) | cart_biomarkers + cart_sequences |
| `scripts/ingest_pubmed.py` | ~5000 | cart_literature |
| `scripts/ingest_clinical_trials.py` | ~1000 | cart_trials |
| `scripts/validate_e2e.py` | -- | End-to-end pipeline validation |
| `scripts/test_rag_pipeline.py` | -- | RAG pipeline integration test |

### 11.4 Seed Data Files

All seed data is stored as JSON files in `data/reference/`:

| File | Records | Description |
|------|---------|-------------|
| assay_seed_data.json | 75 | Cytotoxicity, cytokine, flow, in vivo, TME, exhaustion assay results |
| biomarker_seed_data.json | 60 | Predictive, prognostic, PD, monitoring, resistance biomarkers |
| constructs_seed_data.json | 41 | CAR construct designs, bispecific, armored, universal platforms |
| immunogenicity_biomarker_seed.json | 20 | ADA, HLA, humanization, rejection, cellular immunogenicity biomarkers |
| immunogenicity_sequence_seed.json | 18 | Immunogenicity-related sequence and construct engineering records |
| literature_seed_data.json | 60 | Curated literature records including 2024-2025 landmark papers |
| manufacturing_seed_data.json | 56 | Transduction, expansion, cryo, release testing, POC, iPSC records |
| patent_seed_data.json | 45 | CAR-T patent literature records |
| realworld_seed_data.json | 54 | CIBMTR, institutional, claims-based, disparities, cost RWE studies |
| regulatory_seed_data.json | 40 | FDA/EMA/PMDA/NMPA/Health Canada approval milestones and designations |
| safety_seed_data.json | 71 | CRS, ICANS, cytopenia, infection, secondary malignancy, solid tumor events |
| sequence_seed_data.json | 40 | scFv sequences, VHH nanobodies, binding affinity, structural data |
| trials_seed_data.json | 69 | Clinical trial seed records across hematologic and solid tumors |
| **TOTAL** | **630** | |

---

## 12. API Reference

The FastAPI REST API is defined in `api/main.py` and `api/routes/` and runs on port 8522. There are 13 endpoints total (9 in main.py, 2 in routes/events.py, 1 in routes/meta_agent.py, 2 in routes/reports.py).

### 12.1 Endpoints

#### GET /health

Returns service health with collection count and total vector count.

**Response (200):**
```json
{
  "status": "healthy",
  "collections": 11,
  "total_vectors": 3567622
}
```

**Response (503):** Engine not initialized or Milvus unavailable.

#### GET /collections

Returns all collection names and their record counts.

**Response (200):**
```json
{
  "collections": [
    {"name": "cart_literature", "record_count": 5047},
    {"name": "cart_trials", "record_count": 973},
    ...
  ],
  "total": 11
}
```

#### POST /query

Full RAG query: retrieve evidence from Milvus, augment with knowledge graph, synthesize LLM response.

**Request:**
```json
{
  "question": "Why do CD19 CAR-T therapies fail in relapsed B-ALL?",
  "target_antigen": "CD19",
  "collections": ["cart_literature", "cart_trials"],
  "year_min": 2018,
  "year_max": 2026
}
```

**Response (200):**
```json
{
  "question": "...",
  "answer": "...(LLM-generated response with citations)...",
  "evidence": [
    {
      "collection": "Literature",
      "id": "12345678",
      "score": 0.92,
      "text": "...",
      "metadata": {"title": "...", "year": 2023}
    }
  ],
  "knowledge_context": "## Target Antigen: CD19\n...",
  "collections_searched": 2,
  "search_time_ms": 142.5
}
```

#### POST /search

Evidence-only retrieval (no LLM). Same request schema as /query, returns `SearchResponse` without the `answer` field.

#### POST /find-related

Cross-collection entity linking. Finds all evidence related to an entity across all 10 collections.

**Request:**
```json
{
  "entity": "Yescarta",
  "top_k": 5
}
```

**Response (200):**
```json
{
  "entity": "Yescarta",
  "results": {
    "cart_literature": [...],
    "cart_trials": [...],
    "cart_constructs": [...],
    "cart_safety": [...]
  },
  "total_hits": 18
}
```

#### GET /knowledge/stats

Returns knowledge graph statistics.

**Response (200):**
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

#### GET /metrics

Prometheus-compatible metrics in text exposition format.

**Response (200):**
```
# HELP cart_api_requests_total Total API requests
# TYPE cart_api_requests_total counter
cart_api_requests_total 42

# HELP cart_collection_vectors Number of vectors per collection
# TYPE cart_collection_vectors gauge
cart_collection_vectors{collection="cart_literature"} 5047
...
```

### 12.2 Middleware

- **CORS**: Restricted to 3 origins (`http://localhost:8080`, `http://localhost:8521`, `http://localhost:8522`)
- **Lifespan**: Engine initialization on startup, disconnection on shutdown

### 12.3 Error Handling

- **503**: Engine not initialized, Milvus unavailable, LLM/embedder not loaded
- **500**: Internal query/search/find-related failures (with detail message)
- All errors increment `_metrics["errors_total"]`

---

## 13. UI Guide

The Streamlit UI (defined in `app/cart_ui.py`) runs on port 8521.

### 13.1 Layout

Three main tabs:
1. **Chat**: Conversational RAG interface with streaming responses
2. **Knowledge Graph**: Interactive pyvis graph visualization of entities
3. **Image Analysis**: Upload slide/document images for claim extraction and verification

### 13.2 Sidebar Controls

| Control | Type | Default | Description |
|---------|------|---------|-------------|
| Deep Research Mode | Toggle | Off | Enables autonomous agent with sub-question decomposition |
| Target Antigen Filter | Selectbox | All Targets | Filter search to specific antigen (16 options) |
| Development Stage | Selectbox | All Stages | Filter by CAR-T stage (5 stages) |
| From Year / To Year | Number inputs | 2010-2026 | Temporal date range |
| Apply date filter | Checkbox | Off | Enable/disable date filtering |
| Collection toggles | 11 checkboxes | All on | Enable/disable individual collections (shows live counts) |
| Demo Queries | 13 buttons | -- | Pre-loaded example queries |

### 13.3 Demo Queries

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

### 13.4 Chat Features

- **Streaming responses**: Tokens appear in real-time with a cursor indicator
- **Conversation memory**: Up to 3 prior exchanges injected as context for follow-up queries
- **Comparative detection**: Automatically detects "vs"/"compare" queries and runs dual-entity retrieval
- **Evidence panel**: Expandable panel showing all evidence cards grouped by collection
- **Evidence cards**: Color-coded by collection with relevance badges (high/medium/low), clickable PubMed/ClinicalTrials.gov links
- **Export buttons**: Download as Markdown, JSON, or PDF after each response
- **Deep Research indicator**: Shows search strategy, identified targets, stages, sub-questions, and evidence quality assessment

### 13.5 NVIDIA Theme

Custom CSS implementing the NVIDIA black + green visual identity:
- Background: `#0a0a0f`
- Sidebar: `#12121a`
- NVIDIA Green: `#76B900`
- Cards: `#1a1a24` with `#222230` borders
- Collection badges: 11 unique colors (Literature blue, Trial purple, Construct green, etc.)
- Relevance classes: high (green), medium (yellow), low (gray)

### 13.6 Knowledge Graph Tab

Interactive pyvis network visualization with 5 entity types:
- **Target Antigens**: Green nodes linked to disease (blue), product (purple), and resistance (red) nodes
- **Toxicities**: Red nodes linked to biomarker (teal) and management (indigo) nodes
- **Manufacturing**: Orange nodes linked to parameter (yellow) nodes
- **Biomarkers**: Color-coded by type (predictive/pharmacodynamic/monitoring/resistance)
- **Regulatory**: Indigo nodes linked to indication and designation nodes

Cross-collection entity search is available below the graph.

### 13.7 Image Analysis Tab

Upload PNG/JPG/PDF images for automatic claim extraction:
1. Claude Vision extracts claims and data points from the image
2. Each claim generates a search query
3. Evidence is retrieved across all collections
4. Results show supporting/missing evidence per claim

---

## 14. Metrics & Monitoring

Defined in `src/metrics.py`. Uses `prometheus_client` when available, with no-op stubs as fallback.

### 14.1 Prometheus Metrics

All metrics use the `cart_` prefix.

**Histograms:**
| Metric | Labels | Buckets | Description |
|--------|--------|---------|-------------|
| `cart_query_latency_seconds` | query_type | 0.1, 0.5, 1, 2, 5, 10, 30 | Query processing time |
| `cart_evidence_count` | -- | 0, 5, 10, 15, 20, 25, 30 | Evidence items per query |

**Counters:**
| Metric | Labels | Description |
|--------|--------|-------------|
| `cart_queries_total` | query_type, status | Total queries processed |
| `cart_collection_hits_total` | collection | Hits by collection |
| `cart_llm_tokens_total` | direction | LLM tokens used |

**Gauges:**
| Metric | Labels | Description |
|--------|--------|-------------|
| `cart_active_connections` | -- | Active connections |
| `cart_collection_size` | collection | Records per collection |
| `cart_last_ingest_timestamp` | source | Last ingest timestamp |

### 14.2 Helper Functions

| Function | Purpose |
|----------|---------|
| `record_query(query_type, latency, hit_count, status)` | Record latency + count + status for one query |
| `record_collection_hits(hits_by_collection)` | Increment per-collection hit counters |
| `update_collection_sizes(stats)` | Set current record counts |
| `get_metrics_text()` | Return Prometheus exposition text |

### 14.3 API Metrics Endpoint

The `/metrics` endpoint on the FastAPI server provides:
- Request counters (total, query, search, find-related, errors)
- Per-collection vector counts (when Milvus is available)

### 14.4 Grafana Integration

Metrics are designed to integrate with the HCLS AI Factory Grafana + Prometheus stack:
- Prometheus scrapes on port 9099
- Node Exporter on port 9100
- DCGM Exporter on port 9400
- CAR-T metrics available via `/metrics` on port 8522

---

## 15. Scheduler

Defined in `src/scheduler.py`. Uses APScheduler's `BackgroundScheduler` for automated data refresh.

### 15.1 IngestScheduler

```python
scheduler = IngestScheduler(
    collection_manager=manager,
    embedder=embedder,
    interval_hours=168,  # Weekly (7 days)
)
scheduler.start()   # Non-blocking, runs in daemon thread
scheduler.stop()    # Graceful shutdown
scheduler.get_status()  # Returns next_run_time, last_run_time, job_count
```

### 15.2 Scheduled Jobs

| Job ID | Name | Trigger | Default Interval | Action |
|--------|------|---------|-----------------|--------|
| `refresh_pubmed` | PubMed CAR-T Literature Refresh | Interval | 168 hours (weekly) | Runs PubMedIngestPipeline |
| `refresh_clinical_trials` | ClinicalTrials.gov CAR-T Refresh | Interval | 168 hours (weekly) | Runs ClinicalTrialsIngestPipeline |

### 15.3 Job Execution

Each job:
1. Creates a short-lived ingest pipeline instance
2. Runs the full fetch-parse-embed_and_store pipeline
3. Updates the `cart_last_ingest_timestamp{source=...}` Prometheus gauge
4. Logs completion with record count and elapsed time
5. Catches and logs exceptions without crashing the scheduler

### 15.4 Graceful Degradation

If `apscheduler` is not installed, the module exports a no-op `IngestScheduler` stub that silently ignores all calls, with a warning logged at initialization.

---

## 16. Configuration

All configuration is managed through `CARTSettings` (defined in `config/settings.py`), a Pydantic `BaseSettings` class that reads from environment variables with the `CART_` prefix.

### 16.1 Settings Fields

| Field | Type | Default | Env Var |
|-------|------|---------|---------|
| **Paths** | | | |
| PROJECT_ROOT | Path | (auto-detected) | -- |
| DATA_DIR | Path | `{PROJECT_ROOT}/data` | -- |
| CACHE_DIR | Path | `{DATA_DIR}/cache` | -- |
| REFERENCE_DIR | Path | `{DATA_DIR}/reference` | -- |
| RAG_PIPELINE_ROOT | Path | `/app/rag-chat-pipeline` | `CART_RAG_PIPELINE_ROOT` |
| **Milvus** | | | |
| MILVUS_HOST | str | `localhost` | `CART_MILVUS_HOST` |
| MILVUS_PORT | int | `19530` | `CART_MILVUS_PORT` |
| COLLECTION_LITERATURE | str | `cart_literature` | `CART_COLLECTION_LITERATURE` |
| COLLECTION_TRIALS | str | `cart_trials` | `CART_COLLECTION_TRIALS` |
| COLLECTION_CONSTRUCTS | str | `cart_constructs` | `CART_COLLECTION_CONSTRUCTS` |
| COLLECTION_ASSAYS | str | `cart_assays` | `CART_COLLECTION_ASSAYS` |
| COLLECTION_MANUFACTURING | str | `cart_manufacturing` | `CART_COLLECTION_MANUFACTURING` |
| COLLECTION_GENOMIC | str | `genomic_evidence` | `CART_COLLECTION_GENOMIC` |
| COLLECTION_SAFETY | str | `cart_safety` | `CART_COLLECTION_SAFETY` |
| COLLECTION_BIOMARKERS | str | `cart_biomarkers` | `CART_COLLECTION_BIOMARKERS` |
| COLLECTION_REGULATORY | str | `cart_regulatory` | `CART_COLLECTION_REGULATORY` |
| COLLECTION_SEQUENCES | str | `cart_sequences` | `CART_COLLECTION_SEQUENCES` |
| COLLECTION_REALWORLD | str | `cart_realworld` | `CART_COLLECTION_REALWORLD` |
| **Embeddings** | | | |
| EMBEDDING_MODEL | str | `BAAI/bge-small-en-v1.5` | `CART_EMBEDDING_MODEL` |
| EMBEDDING_DIMENSION | int | `384` | `CART_EMBEDDING_DIMENSION` |
| EMBEDDING_BATCH_SIZE | int | `32` | `CART_EMBEDDING_BATCH_SIZE` |
| **LLM** | | | |
| LLM_PROVIDER | str | `anthropic` | `CART_LLM_PROVIDER` |
| LLM_MODEL | str | `claude-sonnet-4-6` | `CART_LLM_MODEL` |
| ANTHROPIC_API_KEY | str (optional) | None | `CART_ANTHROPIC_API_KEY` |
| **RAG Search** | | | |
| TOP_K_PER_COLLECTION | int | `5` | `CART_TOP_K_PER_COLLECTION` |
| SCORE_THRESHOLD | float | `0.4` | `CART_SCORE_THRESHOLD` |
| WEIGHT_LITERATURE | float | `0.20` | `CART_WEIGHT_LITERATURE` |
| WEIGHT_TRIALS | float | `0.16` | `CART_WEIGHT_TRIALS` |
| WEIGHT_CONSTRUCTS | float | `0.10` | `CART_WEIGHT_CONSTRUCTS` |
| WEIGHT_ASSAYS | float | `0.09` | `CART_WEIGHT_ASSAYS` |
| WEIGHT_MANUFACTURING | float | `0.07` | `CART_WEIGHT_MANUFACTURING` |
| WEIGHT_SAFETY | float | `0.08` | `CART_WEIGHT_SAFETY` |
| WEIGHT_BIOMARKERS | float | `0.08` | `CART_WEIGHT_BIOMARKERS` |
| WEIGHT_REGULATORY | float | `0.06` | `CART_WEIGHT_REGULATORY` |
| WEIGHT_SEQUENCES | float | `0.06` | `CART_WEIGHT_SEQUENCES` |
| WEIGHT_REALWORLD | float | `0.07` | `CART_WEIGHT_REALWORLD` |
| WEIGHT_GENOMIC | float | `0.04` | `CART_WEIGHT_GENOMIC` |
| **PubMed** | | | |
| NCBI_API_KEY | str (optional) | None | `CART_NCBI_API_KEY` |
| PUBMED_MAX_RESULTS | int | `5000` | `CART_PUBMED_MAX_RESULTS` |
| **ClinicalTrials.gov** | | | |
| CT_GOV_BASE_URL | str | `https://clinicaltrials.gov/api/v2` | `CART_CT_GOV_BASE_URL` |
| **API Server** | | | |
| API_HOST | str | `0.0.0.0` | `CART_API_HOST` |
| API_PORT | int | `8522` | `CART_API_PORT` |
| **Streamlit** | | | |
| STREAMLIT_PORT | int | `8521` | `CART_STREAMLIT_PORT` |
| **Prometheus** | | | |
| METRICS_ENABLED | bool | `True` | `CART_METRICS_ENABLED` |
| **Scheduler** | | | |
| INGEST_SCHEDULE_HOURS | int | `168` (weekly) | `CART_INGEST_SCHEDULE_HOURS` |
| INGEST_ENABLED | bool | `False` | `CART_INGEST_ENABLED` |
| **Conversation** | | | |
| MAX_CONVERSATION_CONTEXT | int | `3` | `CART_MAX_CONVERSATION_CONTEXT` |
| **Citation Scoring** | | | |
| CITATION_HIGH_THRESHOLD | float | `0.75` | `CART_CITATION_HIGH_THRESHOLD` |
| CITATION_MEDIUM_THRESHOLD | float | `0.60` | `CART_CITATION_MEDIUM_THRESHOLD` |

### 16.2 Environment Variable Loading

Settings are loaded via:
1. Environment variables with `CART_` prefix (highest priority)
2. `.env` file in the project root (via `env_file=".env"`)
3. Default values (lowest priority)

The Streamlit UI and FastAPI server also independently load `ANTHROPIC_API_KEY` from the rag-chat-pipeline `.env` file as a fallback.

---

## 17. Docker Deployment

### 17.1 Dockerfile (Multi-Stage)

**Stage 1: Builder** (`python:3.10-slim AS builder`)
- Installs build dependencies: build-essential, gcc, g++, libxml2-dev, libxslt1-dev
- Creates virtual environment at `/opt/venv`
- Installs all Python dependencies from requirements.txt

**Stage 2: Runtime** (`python:3.10-slim`)
- Installs minimal runtime libraries: curl, libgomp1, libxml2, libxslt1.1
- Copies virtual environment from builder stage
- Copies application source (config/, src/, app/, scripts/, data/)
- Creates non-root user (`cartuser`)
- Exposes ports 8521 (Streamlit) and 8522 (FastAPI)
- Healthcheck against Streamlit at `/_stcore/health`
- Default CMD: `streamlit run app/cart_ui.py`

### 17.2 Docker Compose Services (6)

| Service | Image | Ports | Purpose |
|---------|-------|-------|---------|
| `milvus-etcd` | quay.io/coreos/etcd:v3.5.5 | -- (internal) | etcd key-value store for Milvus metadata |
| `milvus-minio` | minio/minio:RELEASE.2023-03-20T20-16-18Z | -- (internal) | MinIO object storage for Milvus data |
| `milvus-standalone` | milvusdb/milvus:v2.4-latest | 19530, 9091 | Milvus 2.4 vector database |
| `cart-streamlit` | (built from Dockerfile) | 8521 | Streamlit chat UI |
| `cart-api` | (built from Dockerfile) | 8522 | FastAPI REST server (2 workers) |
| `cart-setup` | (built from Dockerfile) | -- | One-shot: create collections + seed all data |

### 17.3 Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `ANTHROPIC_API_KEY` | Yes | Anthropic API key for Claude |
| `CART_MILVUS_HOST` | No (default: milvus-standalone) | Milvus hostname |
| `CART_MILVUS_PORT` | No (default: 19530) | Milvus port |

### 17.4 Volumes

| Volume | Mount | Purpose |
|--------|-------|---------|
| `etcd_data` | /etcd | etcd persistent storage |
| `minio_data` | /minio_data | MinIO object storage |
| `milvus_data` | /var/lib/milvus | Milvus data directory |

### 17.5 Network

All services connect via the `cart-network` bridge network.

### 17.6 Setup Sequence

The `cart-setup` service runs the following scripts in order:
1. `setup_collections.py --drop-existing --seed-constructs`
2. `seed_knowledge.py`
3. `seed_assays.py`
4. `seed_manufacturing.py`
5. `seed_safety.py`
6. `seed_biomarkers.py`
7. `seed_regulatory.py`
8. `seed_sequences.py`
9. `seed_realworld.py`

After completion, the setup container exits (restart policy: `"no"`).

### 17.7 Deployment Commands

```bash
# Start all services
cp .env.example .env  # Add ANTHROPIC_API_KEY
docker compose up -d

# Watch setup progress
docker compose logs -f cart-setup

# Verify health
curl http://localhost:8522/health
curl http://localhost:8521/_stcore/health

# Access UI
open http://localhost:8521
```

---

## 18. Testing

The test suite contains 415 tests across 7 test files in the `tests/` directory (4,321 lines total).

### 18.1 Test Files

| File | Description |
|------|-------------|
| `tests/conftest.py` | Shared fixtures: mock_embedder, mock_llm_client, mock_collection_manager, sample_search_hits, sample_evidence |
| `tests/test_models.py` | Pydantic model validation, enum values, to_embedding_text(), field constraints |
| `tests/test_knowledge.py` | Knowledge graph: get_target_context, get_toxicity_context, get_manufacturing_context, get_all_context_for_query, entity aliases, comparison context |
| `tests/test_query_expansion.py` | Expansion maps: expand_query, expand_query_by_category, get_expansion_stats, keyword coverage |
| `tests/test_rag_engine.py` | RAG engine: retrieve, query, query_stream, find_related, comparative analysis, citation formatting, prompt building |
| `tests/test_agent.py` | Agent: search_plan, evaluate_evidence, run, generate_report |
| `tests/test_export.py` | Export: export_markdown, export_json, export_pdf, collection-specific table formatting |
| `tests/test_integration.py` | Integration tests: end-to-end pipeline validation, cross-collection queries |

### 18.2 Fixtures

| Fixture | Scope | Description |
|---------|-------|-------------|
| `mock_embedder` | Function | Returns 384-dim zero vectors for any input |
| `mock_llm_client` | Function | Returns "Mock response" for generate(), yields ["Mock ", "response"] for generate_stream() |
| `mock_collection_manager` | Function | Empty search results, dummy stats (42 per collection), no-op connect/disconnect |
| `sample_search_hits` | Function | 5 SearchHit objects across Literature, Trial, Construct, Safety, Manufacturing |
| `sample_evidence` | Function | CrossCollectionResult with 5 hits, knowledge context, and metrics |

### 18.3 Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_models.py -v

# Run with coverage
pytest tests/ --cov=src --cov-report=html
```

### 18.4 Test Coverage Areas

- **Model validation**: Field constraints, enum membership, max_length enforcement, default values, embedding text generation
- **Knowledge graph**: All 3 dictionaries return non-empty context, entity alias resolution, comparison context generation, stats correctness
- **Query expansion**: Keyword matching, category grouping, term deduplication, expansion stats accuracy
- **RAG engine**: End-to-end retrieve flow, parallel search mocking, query expansion integration, knowledge augmentation, citation formatting (PubMed, ClinicalTrials.gov, generic), comparative detection and parsing, prompt assembly
- **Agent**: Search plan extraction, strategy classification (broad/targeted/comparative), evidence quality evaluation, sub-question generation, report structure
- **Export**: Markdown structure, JSON schema correctness, PDF generation (bytes output), collection-specific table columns, filter formatting, comparative export

---

## 19. Project Timeline

### Git History (20 Commits)

| # | Hash | Description | Key Changes |
|---|------|------------|-------------|
| 1 | `3dd5c8a` | Initial scaffold + Week 1 | Project structure, base models, first 5 collections |
| 2 | `cda1093` | .gitignore | Standard Python/data gitignore |
| 3 | `5a6daa8` | Apache 2.0 license | LICENSE file |
| 4 | `fa5b41c` | Real data ingest + validation | PubMed/ClinicalTrials.gov parsers, Pydantic validation |
| 5 | `9fae702` | Full RAG pipeline + UI polish | RAG engine, Streamlit UI, knowledge graph |
| 6 | `c078ac7` | README + setup instructions | Documentation and quickstart |
| 7 | `d344c48` | Seed assays (45 records) | Assay seed data and ingest script |
| 8 | `b93d802` | Seed manufacturing (30 records) | Manufacturing seed data and ingest script |
| 9 | `9af1263` | Design doc with benchmarks | Architecture documentation and performance data |
| 10 | `96ad742` | Fix Claude model ID | Corrected LLM_MODEL to claude-sonnet-4-6 |
| 11 | `3210594` | Citation links + evidence panel | Clickable PubMed/CT.gov links, evidence cards |
| 12 | `7b519c9` | Comparative analysis mode | Entity parsing, dual retrieval, comparison prompts |
| 13 | `052c5ec` | README update | Updated documentation |
| 14 | `2925e13` | Design doc v1.2.0 | Revised architecture documentation |
| 15 | `41bd76f` | Markdown + JSON export | Export system (2 formats) |
| 16 | `9f466a3` | PDF export with NVIDIA theme | reportlab PDF generation with NVIDIA styling |
| 17 | `1ed77e8` | 5 new collections | safety, biomarkers, regulatory, sequences, real-world |
| 18 | `ebc05f2` | v2.0.0 with 21 improvements | Major version with parallel search, settings-driven weights, full knowledge augmentation |
| 19 | `2423fe2` | 3 high-priority gaps | Genomics bridge, patents, immunogenicity |
| 20 | `4e6da4a` | UI theme matching | NVIDIA black + green CSS theme |

### Development Phases

**Phase 1 (Week 1):** Scaffold, first 5 collections (literature, trials, constructs, assays, manufacturing), base ingest pipeline, knowledge graph (targets + toxicities + manufacturing), query expansion (6 maps), RAG engine, Streamlit UI.

**Phase 2 (Week 2):** Real data ingest (PubMed + ClinicalTrials.gov), seed scripts (assays, manufacturing), comparative analysis, citation links, evidence panel, Markdown/JSON/PDF export.

**Phase 3 (Week 3):** 5 new collections (safety, biomarkers, regulatory, sequences, real-world), expanded knowledge graph (biomarkers, regulatory, immunogenicity), expanded query expansion (12 maps), FastAPI REST API, metrics, scheduler.

**Phase 4 (Week 4):** v2.0.0 release with 21 improvements: parallel search, settings-driven configuration, full knowledge augmentation, conversation memory, Deep Research mode, image analysis, genomics bridge, patent data, immunogenicity collections, UI theme matching.

---

## 20. File Inventory

### Complete List of 61 Python Files

#### config/ (1 file, 113 lines)

| File | Lines | Description |
|------|-------|-------------|
| `config/settings.py` | 113 | CARTSettings Pydantic configuration |

#### src/ (28 files, 12,944 lines)

| File | Lines | Description |
|------|-------|-------------|
| `src/__init__.py` | -- | Package init |
| `src/agent.py` | 309 | CARTIntelligenceAgent: plan-search-synthesize |
| `src/collections.py` | 1,004 | CARTCollectionManager: 11 schemas + CRUD |
| `src/export.py` | 1,487 | Markdown, JSON, PDF export with NVIDIA theme |
| `src/knowledge.py` | 2,249 | Knowledge graph: 3 dictionaries, 71 entries, 54 aliases |
| `src/metrics.py` | 404 | Prometheus metrics: counters, histograms, gauges |
| `src/models.py` | 484 | Pydantic models: 13 enums, 10 collections, 4 search/agent |
| `src/query_expansion.py` | 1,592 | 12 expansion maps, 229 keywords, 1,961 terms |
| `src/rag_engine.py` | 754 | CARTRAGEngine: retrieval, synthesis, comparative |
| `src/scheduler.py` | 226 | APScheduler: weekly PubMed/trials refresh |
| `src/utils/__init__.py` | -- | Utils package init |
| `src/utils/pubmed_client.py` | -- | NCBI E-utilities client |
| `src/ingest/__init__.py` | -- | Ingest package init |
| `src/ingest/base.py` | 185 | BaseIngestPipeline: fetch-parse-embed_and_store |
| `src/ingest/literature_parser.py` | -- | PubMed literature ingest |
| `src/ingest/clinical_trials_parser.py` | -- | ClinicalTrials.gov ingest |
| `src/ingest/construct_parser.py` | -- | CAR construct data parser |
| `src/ingest/assay_parser.py` | -- | Assay data parser |
| `src/ingest/manufacturing_parser.py` | -- | Manufacturing data parser |
| `src/ingest/safety_parser.py` | -- | Safety/pharmacovigilance parser |
| `src/ingest/biomarker_parser.py` | -- | Biomarker data parser |
| `src/ingest/regulatory_parser.py` | -- | Regulatory milestone parser |
| `src/ingest/sequence_parser.py` | -- | Molecular sequence parser |
| `src/ingest/realworld_parser.py` | -- | Real-world evidence parser |
| `src/ingest/faers_parser.py` | -- | FDA FAERS safety data parser |
| `src/ingest/dailymed_parser.py` | -- | DailyMed label parser |
| `src/ingest/uniprot_parser.py` | -- | UniProt protein data parser |
| `src/ingest/cibmtr_parser.py` | -- | CIBMTR registry data parser |

#### app/ (1 file, 1,162 lines)

| File | Lines | Description |
|------|-------|-------------|
| `app/cart_ui.py` | 1,162 | Streamlit UI: chat, knowledge graph, image analysis |

#### api/ (5 files, 1,033 lines)

| File | Lines | Description |
|------|-------|-------------|
| `api/__init__.py` | -- | API package init |
| `api/main.py` | 588 | FastAPI REST server: 9 endpoints |
| `api/routes/__init__.py` | -- | Routes package init |
| `api/routes/events.py` | 122 | Pipeline event endpoints (2 endpoints) |
| `api/routes/meta_agent.py` | 142 | Meta-agent /ask endpoint (1 endpoint) |
| `api/routes/reports.py` | 180 | Report generation endpoints (2 endpoints) |

#### scripts/ (15 files, 1,686 lines)

| File | Description |
|------|-------------|
| `scripts/setup_collections.py` | Create all 11 collections |
| `scripts/seed_knowledge.py` | Validate knowledge graph |
| `scripts/seed_assays.py` | Seed 75 assay records |
| `scripts/seed_manufacturing.py` | Seed 56 manufacturing records |
| `scripts/seed_safety.py` | Seed 71 safety records |
| `scripts/seed_biomarkers.py` | Seed 60 biomarker records |
| `scripts/seed_regulatory.py` | Seed 40 regulatory records |
| `scripts/seed_sequences.py` | Seed 40 sequence records |
| `scripts/seed_realworld.py` | Seed 54 real-world records |
| `scripts/seed_patents.py` | Seed 45 patent literature records |
| `scripts/seed_immunogenicity.py` | Seed 38 immunogenicity records (20+18) |
| `scripts/ingest_pubmed.py` | Full PubMed ingest |
| `scripts/ingest_clinical_trials.py` | Full ClinicalTrials.gov ingest |
| `scripts/validate_e2e.py` | End-to-end pipeline validation |
| `scripts/test_rag_pipeline.py` | RAG pipeline integration test |

#### tests/ (7 files, 4,321 lines)

| File | Description |
|------|-------------|
| `tests/__init__.py` | Test package init |
| `tests/conftest.py` | Shared fixtures (5 fixtures) |
| `tests/test_models.py` | Model validation tests |
| `tests/test_knowledge.py` | Knowledge graph tests |
| `tests/test_query_expansion.py` | Query expansion tests |
| `tests/test_rag_engine.py` | RAG engine tests |
| `tests/test_agent.py` | Agent architecture tests |
| `tests/test_export.py` | Export system tests |
| `tests/test_integration.py` | Integration tests |

---

## 21. Dependencies

### requirements.txt

| Package | Version | Category | Purpose |
|---------|---------|----------|---------|
| pydantic | >=2.0 | Core | Data validation and model definitions |
| pydantic-settings | >=2.0 | Core | Environment-based configuration |
| loguru | >=0.7.0 | Core | Structured logging |
| pymilvus | >=2.4.0 | Vector DB | Milvus client library |
| sentence-transformers | >=2.2.0 | Embeddings | BGE-small-en-v1.5 model |
| anthropic | >=0.18.0 | LLM | Claude API client |
| streamlit | >=1.30.0 | Web/API | Chat UI framework |
| fastapi | >=0.109.0 | Web/API | REST API framework |
| uvicorn[standard] | >=0.27.0 | Web/API | ASGI server |
| python-multipart | >=0.0.6 | Web/API | File upload support |
| requests | >=2.31.0 | Ingest | HTTP client |
| lxml | >=5.0.0 | Ingest | PubMed XML parsing |
| biopython | >=1.83 | Ingest | NCBI E-utilities |
| apscheduler | >=3.10.0 | Scheduling | Background job scheduler |
| prometheus-client | >=0.20.0 | Monitoring | Prometheus metrics |
| opentelemetry-api | >=1.29.0 | Observability | OpenTelemetry tracing API |
| opentelemetry-sdk | >=1.29.0 | Observability | OpenTelemetry tracing SDK |
| reportlab | >=4.0.0 | Export | PDF report generation |
| pyvis | >=0.3.0 | Visualization | Knowledge graph rendering |
| numpy | >=1.24.0 | Utilities | Numerical operations |
| tqdm | >=4.65.0 | Utilities | Progress bars |
| python-dotenv | >=1.0.0 | Utilities | .env file loading |

---

## 22. Future Roadmap

### 22.1 VAST AI OS Integration

The CAR-T Intelligence Agent is designed for future integration with the VAST AI OS AgentEngine model:
- `CARTIntelligenceAgent.run()` maps to the AgentEngine entry point
- `search_plan()` maps to the Plan phase
- `rag_engine.retrieve()` maps to the Execute phase
- `evaluate_evidence()` maps to the Reflect phase
- `generate_report()` maps to the Report phase

### 22.2 Additional Agents

The `ai_agent_adds/` directory is designed for multiple domain agents:
- **Imaging Intelligence Agent**: Pathology slide analysis with CAR-T knowledge augmentation
- **Biomarker Intelligence Agent**: Biomarker discovery, biological age modeling, and disease trajectory prediction
- **Clinical Trial Intelligence Agent**: Automated trial landscape monitoring

### 22.3 Data Expansion

- **Expanded PubMed coverage**: Currently ~5,000 records; target 25,000+ with automated weekly refresh
- **Patent database integration**: Full USPTO/EPO CAR-T patent corpus
- **FAERS integration**: Real-time FDA adverse event reporting ingestion
- **EHR data**: De-identified electronic health record integration for outcomes analysis
- **Single-cell RNA-seq**: T-cell phenotype data from public repositories (GEO/ArrayExpress)

### 22.4 Technical Improvements

- **GPU-accelerated embedding**: Move from CPU-based sentence-transformers to GPU inference on DGX Spark
- **Milvus GPU index**: Switch from IVF_FLAT to GPU_IVF_FLAT for faster search at scale
- **Multi-turn reasoning**: Extend agent from 2 sub-questions to full chain-of-thought multi-step retrieval
- **Retrieval evaluation**: Implement retrieval accuracy benchmarks against expert-curated question-answer pairs
- **Fine-tuned embeddings**: Domain-adapted BGE-small for CAR-T terminology

### 22.5 Scaling Targets

| Metric | Current | Target |
|--------|---------|--------|
| Total vectors | 3.57M | 10M+ |
| Literature records | 5,047 | 25,000+ |
| Trial records | 973 | 5,000+ |
| Safety records | 40 | 500+ |
| Collections | 11 | 15+ |
| Query latency (p50) | ~500ms | <200ms |
| Concurrent users | 1 | 10+ |

### 22.6 Hardware Target

**NVIDIA DGX Spark**
- GPU: GB10 (NVIDIA Grace Blackwell)
- Memory: 128GB unified LPDDR5x
- CPU: 20 ARM Cortex-A78AE cores (Grace CPU)
- Interconnect: NVLink-C2C (chip-to-chip)
- Price: $4,699
- OS: Ubuntu-based DGX OS

The entire HCLS AI Factory -- genomics pipeline, RAG/Chat (including this CAR-T Intelligence Agent), and drug discovery pipeline -- is designed to run concurrently on this single workstation, democratizing precision medicine and CAR-T intelligence for institutions that cannot afford traditional HPC infrastructure.

---

*Generated by HCLS AI Factory -- CAR-T Intelligence Agent v2.0.0*
*Author: Adam Jones | License: Apache 2.0 | March 2026*
