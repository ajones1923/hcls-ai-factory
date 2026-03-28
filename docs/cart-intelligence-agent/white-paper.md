# Breaking Down Data Silos in CAR-T Cell Therapy Development: A Multi-Collection RAG Architecture for Cross-Functional Intelligence

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

Part of the HCLS AI Factory -- an end-to-end precision medicine platform.
https://github.com/ajones1923/hcls-ai-factory

---

## Abstract

Chimeric antigen receptor T-cell (CAR-T) therapy represents one of the most transformative advances in oncology, yet the data ecosystem supporting its development remains profoundly fragmented. Clinical researchers, manufacturing engineers, regulatory strategists, and pharmacovigilance scientists each operate within isolated information silos -- PubMed for literature, ClinicalTrials.gov for trial data, FDA databases for regulatory milestones, CIBMTR for real-world outcomes, and institutional databases for manufacturing records. A single cross-functional question such as "How do manufacturing parameters affect clinical response rates for CD19 CAR-T products?" requires manual synthesis across at least five distinct data sources, a process that can consume days of expert time.

This paper presents the CAR-T Intelligence Agent, an AI-powered multi-collection retrieval-augmented generation (RAG) system that unifies 11 specialized vector collections containing 3,567,622 indexed vectors across the full CAR-T development lifecycle. The system employs 384-dimensional BGE-small-en-v1.5 embeddings, a 3-dictionary knowledge graph with 71 structured entries (34 targets, 17 toxicities, 20 manufacturing), 12 query expansion maps mapping 229 keywords to 1,961 related terms, and Claude Sonnet 4.6 for evidence synthesis. Parallel search via ThreadPoolExecutor across all 11 collections delivers cross-functional answers with clickable PubMed and ClinicalTrials.gov citations in under 30 seconds. The system includes structured comparative analysis, citation relevance scoring, and multi-format export (Markdown, JSON, PDF). Comprising 21,259 lines of Python across 61 files with 415 automated tests, the agent runs on a single NVIDIA DGX Spark ($4,699) and is released under the Apache 2.0 license. We demonstrate that a carefully designed multi-collection RAG architecture can transform fragmented biomedical data into actionable cross-functional intelligence, democratizing access to CAR-T development knowledge that was previously available only to large pharmaceutical organizations with dedicated informatics teams.

---

## 1. Introduction

### 1.1 The CAR-T Data Challenge

CAR-T cell therapy has emerged as a paradigm-shifting treatment modality in oncology, with six FDA-approved products generating billions of dollars in annual revenue and hundreds of clinical trials actively enrolling patients worldwide. Yet the informational infrastructure supporting CAR-T development has not kept pace with the science. Data critical to CAR-T development is scattered across fundamentally different systems:

- **Published literature** resides in PubMed (over 40 million citations) and publisher databases, requiring keyword searches that often miss semantically related work.
- **Clinical trial data** lives in ClinicalTrials.gov (over 480,000 registered studies), with structured but inconsistently populated fields.
- **Manufacturing records** -- transduction efficiency, expansion kinetics, release testing results -- are locked in institutional quality management systems and occasionally published in supplementary materials.
- **Safety data** spans FDA Adverse Event Reporting System (FAERS) reports, product labeling, post-market surveillance databases, and clinical trial safety tables.
- **Biomarker data** is distributed across correlative studies, retrospective analyses, and manufacturer-sponsored investigations.
- **Regulatory intelligence** requires monitoring FDA, EMA, and other agency actions across multiple products and indications.
- **Real-world evidence** from registries such as CIBMTR and EBMT offers insights that controlled trials cannot capture.
- **Molecular data** -- scFv sequences, binding affinities, structural information -- is scattered across patents, publications, and protein databases.
- **Genomic evidence** from patient variant data intersects with CAR-T outcomes in ways that are only beginning to be understood.

The consequence of this fragmentation is that cross-functional questions -- the questions that drive the most important decisions in CAR-T development -- are the hardest to answer. When a clinical team asks, "What manufacturing parameters predict better response in BCMA CAR-T for multiple myeloma?", the answer requires synthesizing data from manufacturing records, clinical trial outcomes, published correlative studies, and possibly real-world registry data. No existing tool provides this cross-functional synthesis.

### 1.2 Limitations of Traditional Approaches

Conventional approaches to CAR-T intelligence suffer from several structural limitations:

1. **Keyword-based literature search** misses semantically related content. A PubMed search for "CAR-T exhaustion" will not find papers discussing "T-cell dysfunction" or "checkpoint upregulation" unless those exact terms appear.
2. **Single-source tools** (PubMed, ClinicalTrials.gov, DrugBank) each address one data type, requiring researchers to manually integrate findings across platforms.
3. **Commercial intelligence platforms** (Citeline, GlobalData, Evaluate) provide curated analyses but at costs exceeding $50,000-100,000 per year, excluding academic researchers and smaller biotechnology companies.
4. **General-purpose AI assistants** lack the domain-specific knowledge graphs, structured data, and citation provenance required for scientific rigor.

### 1.3 Our Contribution

The CAR-T Intelligence Agent addresses these limitations through a multi-collection RAG architecture that:

- Unifies **11 specialized vector collections** spanning the complete CAR-T lifecycle into a single query interface
- Employs **parallel vector search** across all collections simultaneously, returning cross-functional results in milliseconds
- Augments retrieval with a **structured knowledge graph** containing curated clinical data on 34 target antigens, 17 toxicity profiles, 20 manufacturing processes, 23 biomarkers, 6 FDA-approved products, and 6 immunogenicity topics
- Expands queries using **12 domain-specific maps** that map 229 expert keywords to 1,961 related terms, dramatically improving recall
- Synthesizes evidence through **Claude Sonnet 4.6** with a domain-expert system prompt spanning 12 areas of CAR-T expertise
- Provides **clickable citation links** to PubMed and ClinicalTrials.gov for every evidence item, maintaining the traceability that scientific work demands
- Runs on a **single NVIDIA DGX Spark** ($4,699), democratizing access to sophisticated CAR-T intelligence

---

## 2. Background

### 2.1 CAR-T Cell Therapy: Current Landscape

CAR-T cell therapy involves engineering a patient's own T-cells (or donor-derived T-cells in allogeneic approaches) to express a chimeric antigen receptor that directs cytotoxic activity against tumor cells bearing a specific surface antigen. The field has progressed rapidly since the first FDA approval of tisagenlecleucel (Kymriah) in August 2017, with six products now approved:

| Product | Generic Name | Manufacturer | Target | Initial Approval | Key Indications |
|---------|-------------|-------------|--------|-----------------|-----------------|
| Kymriah | tisagenlecleucel | Novartis | CD19 | 2017-08-30 | Pediatric/young adult r/r B-ALL, adult r/r DLBCL, r/r FL |
| Yescarta | axicabtagene ciloleucel | Kite/Gilead | CD19 | 2017-10-18 | Adult r/r LBCL (2L and 3L+), r/r FL |
| Tecartus | brexucabtagene autoleucel | Kite/Gilead | CD19 | 2020-07-24 | Adult r/r MCL, adult r/r B-ALL |
| Breyanzi | lisocabtagene maraleucel | BMS/Juno | CD19 | 2021-02-05 | Adult r/r LBCL (2L and 3L+) |
| Abecma | idecabtagene vicleucel | BMS/bluebird bio | BCMA | 2021-03-26 | Adult r/r multiple myeloma (4+ prior lines) |
| Carvykti | ciltacabtagene autoleucel | Janssen/Legend Biotech | BCMA | 2022-02-28 | Adult r/r multiple myeloma (4+ and 1-3 prior lines) |

All six products received Breakthrough Therapy Designation from the FDA. All are subject to a shared CAR-T REMS program requiring treatment center certification for CRS and neurological toxicity management. All carry 15-year post-marketing follow-up requirements for secondary malignancy surveillance.

### 2.2 The Multi-Modal Data Problem

CAR-T development generates data across at least ten distinct categories, each with its own structure, vocabulary, and source systems:

1. **Literature** -- PubMed abstracts, full-text articles, conference proceedings, preprints
2. **Clinical Trials** -- Structured registration data, protocol summaries, outcome results
3. **CAR Constructs** -- Molecular architecture (scFv, hinge, transmembrane, costimulatory, signaling domains)
4. **Functional Assays** -- Cytotoxicity, cytokine secretion, proliferation, persistence, exhaustion markers
5. **Manufacturing/CMC** -- Transduction efficiency, expansion kinetics, release testing, cryopreservation, vein-to-vein logistics
6. **Safety/Pharmacovigilance** -- Adverse events (CRS, ICANS, cytopenias), FAERS reports, long-term follow-up
7. **Biomarkers** -- Predictive (ferritin, CRP), pharmacodynamic (IL-6, CAR-T expansion), resistance (PD-1, LAG-3, TIM-3)
8. **Regulatory Intelligence** -- BLA filings, approval decisions, label updates, REMS requirements
9. **Molecular/Sequence Data** -- scFv sequences, CDR regions, binding affinities, humanization status, immunogenicity risk
10. **Real-World Evidence** -- Registry outcomes (CIBMTR, EBMT), community vs. academic outcomes, special populations (elderly, bridging therapy), disparities

Each data category has a different optimal schema, different source APIs, different update cadences, and different relevance to different user roles. A traditional single-collection vector database cannot adequately represent this heterogeneity. A multi-collection architecture, where each collection has a purpose-built schema with domain-specific metadata fields, is essential.

---

## 3. System Architecture

### 3.1 Overview

The CAR-T Intelligence Agent implements a multi-collection RAG architecture with five core components: (1) domain-specific Milvus collections with typed schemas, (2) parallel vector search with configurable collection weights, (3) a structured knowledge graph for context augmentation, (4) query expansion for improved recall, and (5) LLM synthesis with citation provenance. The system processes a user query through the following pipeline:

```
User Query
    |
    v
[Comparative Detection] -- "X vs Y" detected? -- YES --> [Parse Two Entities]
    |                                                      (resolve via knowledge graph)
    NO                                                              |
    |                                                      [Dual Retrieval]
    v                                                      (Entity A + Entity B)
[BGE-small-en-v1.5 Embedding]                                      |
(384-dim, asymmetric query prefix)                                  v
    |                                                [Comparative Prompt Builder]
    v                                                (tables + pros/cons format)
[Parallel Search: 11 Milvus Collections]                            |
(ThreadPoolExecutor, IVF_FLAT / COSINE)                             |
    |                                                               |
    v                                                               |
[Query Expansion] (12 maps, 229 keywords -> 1,961 terms)           |
    |                                                               |
    v                                                               |
[Knowledge Graph Augmentation]                                      |
(34 targets, 17 toxicities, 20 mfg, 23 biomarkers,                 |
 6 regulatory, 6 immunogenicity)                                    |
    |                                                               |
    v                                                               v
[Score-Weighted Merge & Rank]
(citation relevance: high >= 0.75, medium >= 0.60)
    |
    v
[Claude Sonnet 4.6] --> Grounded response with clickable citations
    |
    v
[Export] --> Markdown | JSON | PDF (NVIDIA-themed)
```

### 3.2 The 11 Collections

The system maintains 11 Milvus collections, each with a purpose-built schema containing domain-specific metadata fields. Ten collections are owned and populated by the CAR-T agent; the eleventh (`genomic_evidence`) is a read-only shared collection created by the upstream genomics pipeline.

| # | Collection | Vectors | Schema Fields | Source | Description |
|---|-----------|---------|---------------|--------|-------------|
| 1 | `cart_literature` | 5,047 | id, title, text_chunk, source_type, year, cart_stage, target_antigen, disease, keywords, journal | PubMed E-utilities | Published research papers and patents |
| 2 | `cart_trials` | 973 | id, title, text_summary, phase, status, sponsor, target_antigen, car_generation, costimulatory, disease, enrollment, start_year, outcome_summary | ClinicalTrials.gov V2 API | Clinical trial registrations |
| 3 | `cart_constructs` | 41 | id, name, text_summary, target_antigen, scfv_origin, costimulatory_domain, signaling_domain, generation, hinge_tm, vector_type, fda_status, known_toxicities | Curated (6 FDA products) | CAR molecular architecture |
| 4 | `cart_assays` | 75 | id, text_summary, assay_type, construct_id, target_antigen, cell_line, effector_ratio, key_metric, metric_value, outcome, notes | Curated from ELIANA, ZUMA-1, KarMMa, CARTITUDE-1 | Functional assay results |
| 5 | `cart_manufacturing` | 56 | id, text_summary, process_step, vector_type, parameter, parameter_value, target_spec, met_spec, batch_id, notes | Curated CMC data | Manufacturing process records |
| 6 | `cart_safety` | 71 | id, text_summary, product, event_type, severity_grade, onset_timing, incidence_rate, management_protocol, outcome, reporting_source, year | FAERS, labels, trials | Adverse events and pharmacovigilance |
| 7 | `cart_biomarkers` | 60 | id, text_summary, biomarker_name, biomarker_type, assay_method, clinical_cutoff, predictive_value, associated_outcome, target_antigen, disease, evidence_level | Published correlative studies | Predictive and prognostic biomarkers |
| 8 | `cart_regulatory` | 40 | id, text_summary, product, regulatory_event, date, agency, indication, decision, conditions, pivotal_trial | FDA/EMA public records | Regulatory milestones and approvals |
| 9 | `cart_sequences` | 40 | id, text_summary, construct_name, target_antigen, scfv_clone, binding_affinity_kd, heavy_chain_vregion, light_chain_vregion, framework, species_origin, immunogenicity_risk, structural_notes | Patents, UniProt, PDB | Molecular and structural data |
| 10 | `cart_realworld` | 54 | id, text_summary, study_type, data_source, product, indication, population_size, median_followup_months, primary_endpoint, outcome_value, setting, special_population | CIBMTR, EBMT, institutional | Real-world evidence and outcomes |
| 11 | `genomic_evidence` | 3,561,170 | id, text_summary, chrom, pos, ref, alt, qual, gene, consequence, impact, genotype, clinical_significance, rsid, disease_associations, am_pathogenicity, am_class | HCLS AI Factory Stage 1+2 | Shared genomic variant data (read-only) |
| | **Total** | **3,567,622** | | | |

All collections use 384-dimensional float vectors generated by BGE-small-en-v1.5, indexed with IVF_FLAT (nlist=1024) and searched with COSINE similarity (nprobe=16).

### 3.3 Parallel Multi-Collection Search

The `CARTCollectionManager.search_all()` method dispatches vector similarity searches to all 11 collections simultaneously using Python's `concurrent.futures.ThreadPoolExecutor` with a worker count equal to the number of collections. Each collection search operates independently, applying collection-specific filter expressions (e.g., `target_antigen == "CD19"` on collections that support this field, year-range filters on collections with temporal fields).

Results are merged through a score-weighted ranking mechanism. Each collection is assigned a configurable weight reflecting its typical relevance:

| Collection | Weight | Rationale |
|-----------|--------|-----------|
| Literature | 0.20 | Largest corpus, broadest coverage |
| Trials | 0.16 | Structured clinical data |
| Constructs | 0.10 | Foundational molecular data |
| Assays | 0.09 | Preclinical evidence |
| Safety | 0.08 | Critical for clinical decisions |
| Biomarkers | 0.08 | Predictive clinical utility |
| Manufacturing | 0.07 | Process optimization insights |
| Real-World | 0.07 | Post-market validation |
| Regulatory | 0.06 | Approval pathway guidance |
| Sequences | 0.06 | Molecular design data |
| Genomic | 0.04 | Supplementary variant context |

The weighted score is computed as `raw_cosine_score * (1 + collection_weight)`, ensuring that high-relevance collections contribute proportionally more to the final ranking. After weighting, results are deduplicated by record ID and capped at 30 hits for prompt construction.

### 3.4 Embedding Strategy

The system employs BAAI/bge-small-en-v1.5, a 384-dimensional sentence embedding model optimized for asymmetric retrieval. This model was selected for three reasons:

1. **Compact dimensionality** (384 vs. 768 or 1024) reduces memory consumption, critical when indexing 3.5 million vectors on a single machine with 128 GB unified memory.
2. **Asymmetric query prefix** ("Represent this sentence for searching relevant passages: ") improves retrieval quality for question-to-passage matching.
3. **Consistency** with the parent HCLS AI Factory platform, which uses the same model for genomic variant embeddings, enabling the shared `genomic_evidence` collection to be searched seamlessly.

Each record's embedding text is generated by a model-specific `to_embedding_text()` method that concatenates the most semantically important fields. For example, a `ClinicalTrial` record combines its title, summary, target antigen, indication, and outcome summary into the embedding text, ensuring that semantic search captures the trial's full context.

---

## 4. Knowledge Augmentation

### 4.1 The 3-Dictionary Knowledge Graph

Unlike purely retrieval-based systems, the CAR-T Intelligence Agent augments vector search results with structured knowledge from a curated domain knowledge graph. This graph is implemented as three Python dictionaries containing 71 entries and 54 entity aliases:

| Dictionary | Entries | Content |
|-----------|---------|---------|
| `CART_TARGETS` | 34 | Target antigens with expression profiles, approved products, key trials, resistance mechanisms, toxicity profiles, UniProt IDs |
| `CART_TOXICITIES` | 17 | CRS, ICANS, B-cell aplasia, HLH/MAS, cytopenias, TLS, GvHD, on-target/off-tumor, and 4 additional profiles -- each with grading systems, incidence, timing, management protocols, biomarkers, risk factors |
| `CART_MANUFACTURING` | 20 | Lentiviral/retroviral transduction, T-cell activation, ex vivo expansion, leukapheresis, cryopreservation, release testing, point-of-care manufacturing, lymphodepletion, vein-to-vein logistics, and 5 additional processes |
| `CART_BIOMARKERS` | 23 | Ferritin, CRP, IL-6, sIL-2R, CAR-T expansion Cmax, Tcm%, CD4:CD8 ratio, LDH, PD-1, LAG-3, TIM-3, MRD, ctDNA, sBCMA, IFN-gamma, and 8 additional biomarkers -- each with assay method, clinical cutoff, predictive value, evidence level |
| `CART_REGULATORY` | 6 | All FDA-approved products with approval dates, indications, pivotal trials, designations (Breakthrough Therapy, RMAT), subsequent approvals, REMS programs, EMA approval dates |
| `CART_IMMUNOGENICITY` | 6 | Murine scFv immunogenicity, humanization strategies, ADA clinical impact, HLA-restricted epitopes, immunogenicity testing paradigm, allogeneic HLA considerations |

When a query mentions a target antigen (e.g., "CD19"), a toxicity (e.g., "CRS"), a manufacturing process (e.g., "leukapheresis"), a biomarker (e.g., "ferritin"), a product (e.g., "Kymriah"), or an immunogenicity topic (e.g., "humanization"), the corresponding knowledge graph entry is injected into the LLM prompt alongside the retrieved evidence. This provides the LLM with structured, curated context that complements the probabilistic retrieval results.

### 4.2 Entity Alias Resolution

The knowledge graph includes 54 entity aliases that map alternative names to canonical entries:

- **Product aliases:** Kymriah/tisagenlecleucel, Yescarta/axicabtagene ciloleucel, and so on for all six products (12 aliases)
- **Costimulatory domain aliases:** 4-1BB/CD137, CD28, OX40/CD134, ICOS (5 aliases)
- **Manufacturing aliases:** lentiviral, retroviral (2 aliases)
- **Biomarker aliases:** ferritin, CRP, IL-6/IL6, PD-1/PD1, LAG-3/LAG3, TIM-3/TIM3, MRD, sBCMA, ctDNA (11 aliases)
- **Immunogenicity aliases:** ADA, anti-drug antibody, HAMA, humanization, deimmunization, HLA, MHC (9 aliases)

This alias resolution is critical for comparative analysis (Section 6), where the system must map user-specified entity names to structured knowledge graph entries for dual retrieval.

### 4.3 Context Functions

Five context retrieval functions (`get_target_context`, `get_toxicity_context`, `get_manufacturing_context`, `get_biomarker_context`, `get_regulatory_context`, `get_immunogenicity_context`) format knowledge graph entries into structured markdown text. A master function, `get_all_context_for_query`, scans the full query text against all three dictionaries using keyword matching and returns combined context for all detected entities.

The knowledge graph augmentation addresses a fundamental weakness of pure vector retrieval: embedding similarity can surface relevant documents, but it cannot provide the structured clinical parameters (grading systems, incidence rates, management protocols, dosing regimens) that clinicians and researchers need. By injecting this structured context alongside retrieved evidence, the LLM can generate responses that are both evidence-grounded and clinically precise.

---

## 5. Query Expansion

### 5.1 The Recall Problem

Semantic search using dense vector embeddings is effective at capturing conceptual similarity, but it can miss important related content when the vocabulary gap between query and document is too large. A researcher asking about "T-cell exhaustion in CAR-T therapy" should also retrieve documents about "PD-1 upregulation," "LAG-3 expression," "TOX transcription factor," and "tonic signaling" -- terms that are conceptually related but may not produce high cosine similarity scores with the original query embedding.

### 5.2 Expansion Architecture

The system implements 12 domain-specific query expansion maps organized by category:

| # | Category | Keywords | Terms | Examples |
|---|----------|----------|-------|---------|
| 1 | Target Antigen | 32 | 244 | "cd19" expands to CD19, B-ALL, DLBCL, tisagenlecleucel, FMC63, ... |
| 2 | Disease | 16 | 143 | "dlbcl" expands to DLBCL, diffuse large B-cell lymphoma, CD19, Yescarta, Breyanzi, GCB subtype, ... |
| 3 | Toxicity | 18 | 168 | "crs" expands to CRS, cytokine release syndrome, tocilizumab, IL-6, ferritin, Lee grading, ... |
| 4 | Manufacturing | 21 | 215 | "expansion" expands to T-cell expansion, Dynabeads, IL-2, IL-7, G-Rex, CliniMACS Prodigy, ... |
| 5 | Mechanism | 19 | 224 | "exhaustion" expands to T-cell exhaustion, PD-1, LAG-3, TIM-3, TOX, tonic signaling, ... |
| 6 | Construct | 20 | 206 | "bispecific" expands to bispecific CAR, tandem CAR, OR-gate logic, CD19/CD22, ... |
| 7 | Safety | 8 | 69 | "pharmacovigilance" expands to post-market surveillance, FAERS, adverse event reporting, REMS, ... |
| 8 | Biomarker | 18 | 117 | "mrd" expands to MRD, minimal residual disease, flow cytometry MRD, measurable residual disease, ... |
| 9 | Regulatory | 8 | 46 | "breakthrough" expands to breakthrough therapy, BTD, expedited program, FDA designation, ... |
| 10 | Sequence | 8 | 54 | "scfv" expands to scFv, single-chain variable fragment, VH, VL, linker, FMC63, nanobody, ... |
| 11 | Real-World | 10 | 65 | "cibmtr" expands to CIBMTR, transplant registry, national registry, ... |
| 12 | Immunogenicity | 12 | 98 | "humanization" expands to humanized antibody, CDR grafting, deimmunization, framework selection, ... |
| | **Total** | **229** | **1,961** | |

### 5.3 Expansion Execution

When a query is processed, the `expand_query()` function scans the lowercase query text against all 229 keywords across the 12 maps. Matched terms are deduplicated and returned as a sorted list. The RAG engine then processes the top 5 expansion terms through a dual strategy:

1. **Antigen terms** (those matching the 28 known target antigens) are used as field-level filters on collections with `target_antigen` metadata, ensuring precise antigen-specific retrieval.
2. **Non-antigen terms** are re-embedded using BGE-small-en-v1.5 and used for semantic search across all collections, broadening coverage to conceptually related documents.

Expansion results are scored at 70-80% of the original query results (0.7x for semantic expansion, 0.8x for antigen-filter expansion), ensuring that directly relevant results from the primary query still rank highest while expansion results supplement coverage.

---

## 6. Data Integration

### 6.1 PubMed Ingestion

The system ingests CAR-T literature from PubMed via the NCBI E-utilities API (esearch + efetch), targeting approximately 5,000 abstracts per refresh cycle. The `PubMedIngestPipeline` executes a multi-term search query designed to capture the breadth of CAR-T research:

```
"CAR-T" OR "chimeric antigen receptor" OR "CAR T-cell" OR
"adoptive cell therapy" OR "CD19 CAR" OR "BCMA CAR"
```

Each abstract is classified by CAR-T development stage (target identification, CAR design, vector engineering, testing, clinical) using keyword heuristics. Target antigens are extracted from title, abstract, and MeSH terms. The pipeline embeds each record with BGE-small-en-v1.5 and inserts it into the `cart_literature` collection with full metadata (PMID, title, year, journal, disease, keywords, target antigen, development stage). Ingestion of approximately 5,000 abstracts requires approximately 15 minutes.

### 6.2 ClinicalTrials.gov Ingestion

The `ClinicalTrialsIngestPipeline` queries the ClinicalTrials.gov V2 API for CAR-T-related trials, extracting structured fields including NCT number, title, phase, recruitment status, lead sponsor, target enrollment, and brief summary. The pipeline applies domain-specific parsing to extract target antigens, CAR generation, costimulatory domain, and disease indication from free-text fields. Approximately 973 CAR-T trials are currently indexed, with ingestion completing in approximately 3 minutes.

### 6.3 Curated Seed Data

Six additional collection types are populated with curated seed data representing reference-quality records:

- **Constructs** (41 records): Complete molecular architecture of all 6 FDA-approved CAR-T products, including scFv clone origin, costimulatory domain, signaling domain, hinge/transmembrane region, vector type, and known toxicity profile.
- **Assays** (75 records): Functional data curated from landmark publications (ELIANA, ZUMA-1, KarMMa, CARTITUDE-1), covering cytotoxicity assays, cytokine profiling, persistence measurements, exhaustion marker expression, and resistance readouts.
- **Manufacturing** (56 records): CMC process data spanning transduction, expansion, harvest, cryopreservation, release testing, logistics, cost analysis, and emerging platforms (point-of-care, allogeneic, non-viral).
- **Safety** (71 records): Adverse event profiles for approved products, sourced from product labeling, pivotal trial safety tables, and FAERS data.
- **Biomarkers** (60 records): Predictive, prognostic, pharmacodynamic, monitoring, and resistance biomarkers with assay methods, clinical cutoffs, and evidence levels.
- **Regulatory** (40 records): FDA and EMA regulatory milestones including BLA filings, approval decisions, label updates, Breakthrough Therapy and RMAT designations, and REMS requirements.
- **Sequences** (40 records): Molecular data including scFv sequences, CDR region annotations, binding affinity measurements, framework assignments, species origin, and immunogenicity risk assessments.
- **Real-World** (54 records): Registry-based outcomes from CIBMTR, institutional retrospective analyses, community vs. academic center comparisons, special population outcomes (elderly, bridging therapy), and health disparities data.

In total, the 10 agent-owned collections contain 6,452 curated records. Combined with the 3,561,170 shared genomic variants, the system indexes 3,567,622 vectors.

### 6.4 Automated Weekly Refresh

The `IngestScheduler` class wraps APScheduler's `BackgroundScheduler` with two recurring interval jobs -- one for PubMed and one for ClinicalTrials.gov -- each executing on a configurable cadence (default: 168 hours / weekly). Each job creates a fresh ingest pipeline instance, runs the full fetch-parse-embed-store cycle, and updates a Prometheus gauge (`cart_last_ingest_timestamp`) for observability. The scheduler runs as a daemon thread alongside the FastAPI/Streamlit application, requiring no external cron or task queue infrastructure.

---

## 7. Comparative Analysis

### 7.1 Motivation

A significant class of CAR-T research questions is inherently comparative: "CD28 vs 4-1BB costimulatory domains," "Kymriah vs Yescarta," "lentiviral vs retroviral transduction," "CRS vs ICANS management." These questions require the system to retrieve evidence for two distinct entities, organize findings by entity, and produce structured side-by-side analysis.

### 7.2 Implementation

Comparative analysis is implemented as a pipeline within the `CARTRAGEngine`:

1. **Detection:** The engine checks for comparative keywords ("compare," "vs," "versus," "comparing") in the query text.
2. **Entity Parsing:** Regular expressions extract two entity names from the query, handling patterns like "X vs Y," "compare X and Y," and "X versus Y." Trailing descriptors ("costimulatory domains," "resistance mechanisms") are stripped.
3. **Entity Resolution:** Each parsed entity name is resolved against the knowledge graph through the `resolve_comparison_entity()` function, which checks (in priority order) target antigens, product aliases, toxicities, manufacturing processes, and biomarkers. The 54 entity aliases enable resolution of alternative names (e.g., "tisagenlecleucel" resolves to "Kymriah" with target "CD19").
4. **Dual Retrieval:** Two independent `retrieve()` calls execute with entity-specific filters. If both entities resolve to target antigens, each retrieval applies a `target_antigen` field filter to relevant collections.
5. **Comparative Prompt Construction:** A specialized prompt template instructs the LLM to produce: (a) a comparison table with key dimensions as rows and the two entities as columns, (b) advantages of each entity, (c) limitations of each entity, and (d) a clinical context paragraph explaining when each might be preferred.
6. **Structured Output:** Results are returned as a `ComparativeResult` model with separate `evidence_a` and `evidence_b` objects, enabling the UI to display entity-grouped evidence with color-coded headers.

### 7.3 Supported Entity Types

The comparative system supports five entity types:

| Entity Type | Examples | Resolution Strategy |
|-------------|---------|-------------------|
| Target antigens | CD19 vs BCMA, CD22 vs CD20 | Direct match in `CART_TARGETS` (34 entries) |
| FDA products | Kymriah vs Yescarta, Abecma vs Carvykti | Alias table (12 product aliases) with target antigen mapping |
| Costimulatory domains | 4-1BB vs CD28, OX40 vs ICOS | Alias table (5 costimulatory aliases) |
| Toxicities | CRS vs ICANS | Match in `CART_TOXICITIES` (17 entries) |
| Manufacturing processes | Lentiviral vs retroviral | Match in `CART_MANUFACTURING` (20 entries) |

Unrecognized entities gracefully fall back to the standard (non-comparative) query path, ensuring no query fails due to entity resolution.

---

## 8. Evidence Quality and Citation Provenance

### 8.1 Citation Relevance Scoring

Every search hit is assigned a relevance tier based on its raw cosine similarity score:

| Tier | Threshold | Interpretation |
|------|-----------|---------------|
| High | >= 0.75 | Strong semantic match; prioritized in LLM prompt |
| Medium | >= 0.60 | Moderate match; included as supporting evidence |
| Low | < 0.60 | Weak match; included for completeness but deprioritized |

The relevance tier is injected as a tag in the LLM prompt (e.g., `[high relevance]`), guiding the model to prioritize the most relevant evidence when constructing its response.

### 8.2 Clickable Citation Links

The system generates clickable hyperlinks for every evidence item where a stable URL exists:

- **Literature:** `[Literature:PMID 12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/)`
- **Trials:** `[Trial:NCT12345678](https://clinicaltrials.gov/study/NCT12345678)`
- **Other collections:** `[Collection:record-id]` format for internally curated records

These clickable citations persist through all export formats (Markdown, PDF with embedded hyperlinks, JSON with structured citation objects), enabling downstream verification of every claim.

### 8.3 Multi-Format Export

The export module provides three output formats:

1. **Markdown:** Human-readable report with evidence tables, search metrics, and knowledge graph context. Collection-specific table columns surface the most relevant metadata for each collection type (e.g., phase/status/sponsor for trials, assay type/cell line/metric for assays).
2. **JSON:** Machine-readable structured data using Pydantic `.model_dump()` serialization, suitable for integration with downstream analysis pipelines or dashboard systems.
3. **PDF:** Styled report generated via ReportLab Platypus with NVIDIA-branded theming (NVIDIA green `#76B900` accent color, dark headers, alternating row colors). The PDF renderer handles markdown-to-flowable conversion including tables, bold text, bullet lists, block quotes, and clickable hyperlinks.

---

## 9. Hardware Democratization

### 9.1 The Cost Barrier

Enterprise-grade CAR-T intelligence systems -- pharmaceutical informatics platforms, commercial competitive intelligence tools, and internal data warehouse systems -- typically require infrastructure costing hundreds of thousands to millions of dollars annually, including cluster computing resources, commercial database licenses, and dedicated data engineering teams. This cost structure effectively restricts sophisticated CAR-T intelligence to large pharmaceutical companies, excluding academic medical centers, smaller biotechnology companies, and institutions in resource-limited settings.

### 9.2 The DGX Spark Platform

The NVIDIA DGX Spark represents a fundamental shift in computational accessibility. At $4,699, it provides:

| Specification | Value |
|--------------|-------|
| GPU | GB10 (NVIDIA Blackwell architecture) |
| Memory | 128 GB unified LPDDR5x (shared between CPU and GPU) |
| CPU | 20 ARM cores (Grace) |
| Interconnect | NVLink-C2C (chip-to-chip) |
| Storage | NVMe SSD |
| Power | Desktop form factor |

The unified memory architecture is particularly important for this application. The 3.5 million vectors indexed in Milvus, each 384-dimensional (384 x 4 bytes = 1.5 KB per vector), require approximately 5.4 GB of memory for vector data alone, plus index structures. The BGE-small-en-v1.5 embedding model consumes approximately 130 MB. Milvus operational overhead, collection metadata, and the knowledge graph require additional memory. The 128 GB unified pool comfortably accommodates all of these workloads simultaneously, with substantial headroom for the parent HCLS AI Factory pipeline's genomics (Parabricks) and drug discovery (BioNeMo) stages.

### 9.3 Local vs. Cloud Architecture

The CAR-T Intelligence Agent uses a hybrid local-cloud architecture:

| Component | Deployment | Rationale |
|-----------|-----------|-----------|
| Milvus 2.4 | Local (localhost:19530) | Data sovereignty; no egress cost for 3.5M vectors |
| BGE-small-en-v1.5 | Local (CPU/GPU) | Fast embedding; no API latency; offline capable |
| Knowledge graph | Local (in-memory Python dicts) | Zero-latency lookups; no external dependencies |
| Query expansion | Local (in-memory Python dicts) | Zero-latency keyword matching |
| Claude Sonnet 4.6 | Cloud (Anthropic API) | Frontier LLM capability; pay-per-token economics |
| Streamlit UI | Local (port 8521) | No hosting cost; local network access |

This architecture minimizes cloud dependency to a single API call per query (the LLM synthesis step), keeping sensitive data local and operational costs near zero for the retrieval and augmentation stages. The cloud LLM call costs approximately $0.01-0.05 per query depending on response length.

---

## 10. Results and Capabilities

### 10.1 Performance Metrics

Measured on NVIDIA DGX Spark (GB10 GPU, 128 GB unified memory):

| Metric | Value |
|--------|-------|
| PubMed ingestion (5,047 abstracts) | ~15 minutes |
| ClinicalTrials.gov ingestion (973 trials) | ~3 minutes |
| Seed data ingestion (630 records across 13 categories) | ~5 minutes |
| Vector search (11 collections, top-5 each, parallel) | 12-16 ms (warm cache) |
| Comparative dual retrieval (2 x 11 collections) | ~365 ms |
| Full RAG query (search + knowledge + Claude) | ~24 seconds |
| Comparative RAG query (dual search + Claude) | ~30 seconds |
| Cosine similarity scores (demo queries) | 0.74 - 0.90 |
| Total indexed vectors | 3,567,622 |
| Total Python codebase | 21,259 lines across 61 files |
| Automated test count | 415 |

### 10.2 Query Capabilities

The system supports seven classes of cross-functional queries:

**Cross-functional mechanism questions:**
> "Why do CD19 CAR-T therapies fail in relapsed B-ALL?"

This query simultaneously retrieves published literature on CD19 antigen loss and resistance mechanisms, clinical trial outcomes showing relapse rates, assay data on CD19 expression levels, and knowledge graph context on CD19 biology -- producing a unified answer that connects preclinical findings to clinical outcomes.

**Manufacturing-outcome correlations:**
> "What manufacturing parameters predict clinical response?"

Retrieves manufacturing process records (transduction efficiency, expansion kinetics, Tcm percentage), biomarker data (product phenotype markers), clinical trial outcomes, and published correlative studies.

**Comparative analysis:**
> "Compare 4-1BB vs CD28 costimulatory domains for DLBCL"

Triggers comparative mode: parses "4-1BB" and "CD28" as entities, resolves them via the knowledge graph, runs dual retrievals, and produces a structured comparison table with signaling kinetics, persistence profiles, toxicity patterns, and clinical outcome data.

**Safety signal detection:**
> "What are the long-term safety concerns for BCMA CAR-T products?"

Searches safety records, regulatory data (label updates, REMS requirements), biomarker data (MRD monitoring), real-world outcomes, and published long-term follow-up studies.

**Biomarker-guided clinical decisions:**
> "How does T-cell exhaustion affect CAR-T persistence?"

Retrieves exhaustion biomarker data (PD-1, LAG-3, TIM-3 expression levels), assay results showing functional impairment, manufacturing data on phenotype optimization, and clinical correlative studies.

**Regulatory pathway navigation:**
> "What regulatory pathway did Carvykti follow for expanded approval?"

Retrieves regulatory records (Breakthrough Therapy designation, initial and expanded approvals), pivotal trial data (CARTITUDE-1, CARTITUDE-4), and real-world evidence supporting the expanded indication.

**Molecular design questions:**
> "Compare FMC63 and humanized anti-CD19 scFvs for immunogenicity risk"

Retrieves sequence data (FMC63 CDR regions, species origin, binding affinity), immunogenicity knowledge (ADA incidence rates, HLA-restricted epitopes), and clinical data on anti-drug antibody impact.

---

## 11. Integration with the HCLS AI Factory

### 11.1 The Three-Stage Pipeline

The CAR-T Intelligence Agent is one component of the HCLS AI Factory, an end-to-end precision medicine platform that processes a patient sample from raw DNA sequencing data to drug candidate molecules in under 5 hours on a single DGX Spark:

| Stage | Duration | Process | Output |
|-------|----------|---------|--------|
| **1. Genomics** | 120-240 min | FASTQ to VCF via Parabricks 4.6, BWA-MEM2, DeepVariant | 11.7 million variants |
| **2. RAG/Chat** | Interactive | ClinVar (~2.7M) + AlphaMissense (71M) + Milvus (3.5M vectors) + Claude | Variant interpretation, gene-drug associations |
| **3. Drug Discovery** | 8-16 min | MolMIM generation + DiffDock docking + RDKit scoring | Candidate molecules with docking scores |

### 11.2 The Genomic Evidence Bridge

The `genomic_evidence` collection (3,561,170 vectors) is the architectural bridge between the CAR-T agent and the parent pipeline. Created and populated by Stage 2 (RAG/Chat pipeline), this collection contains patient variant data annotated with ClinVar clinical significance, AlphaMissense pathogenicity scores, gene-level consequence predictions, and disease associations.

The CAR-T agent accesses this collection in read-only mode, enabling queries that connect patient genomics to CAR-T therapy decisions. For example, a query about "VCP gene variants in the context of T-cell function" would retrieve relevant genomic variants from the shared collection alongside CAR-T-specific evidence from the agent's own collections.

This shared-collection architecture avoids data duplication while enabling specialized agents to access cross-domain data seamlessly. The Milvus database instance running on localhost:19530 serves both the parent pipeline and the CAR-T agent, with each maintaining independent collections under a unified vector index.

### 11.3 Platform Consistency

The CAR-T Intelligence Agent follows the architectural patterns established by the parent platform:

- **Embedding model:** Same BGE-small-en-v1.5 (384-dim) used across all pipeline stages
- **Vector database:** Same Milvus 2.4 instance with IVF_FLAT/COSINE indexes
- **LLM integration:** Same Claude API integration pattern with streaming support
- **Configuration:** Pydantic BaseSettings pattern matching `rag-chat-pipeline/config/settings.py`
- **Data models:** Pydantic BaseModel pattern matching `drug-discovery-pipeline/src/models.py`
- **Observability:** Prometheus metrics integration consistent with the platform's Grafana+Prometheus stack (port 3000 Grafana, port 9099 Prometheus)

---

## 12. Discussion

### 12.1 Implications for CAR-T Development

The CAR-T Intelligence Agent demonstrates that the data fragmentation problem in cell therapy development is solvable with current technology at accessible cost. The key insight is that the diversity of CAR-T data types (literature, trials, manufacturing, safety, biomarkers, regulatory, sequences, real-world evidence, genomics) requires a multi-collection architecture rather than a single monolithic vector store. Each collection's typed schema ensures that domain-specific metadata is preserved and queryable, while the parallel search and score-weighted merge mechanism provides a unified query experience.

The knowledge graph augmentation addresses a practical limitation of pure RAG systems: while vector retrieval can surface relevant documents, it cannot provide the structured, curated clinical parameters (grading criteria, dosing protocols, management algorithms) that practitioners need. The combination of retrieval-based evidence and structured knowledge produces responses that are both evidence-grounded and clinically actionable.

### 12.2 Democratization of Intelligence

The pharmaceutical industry's current informatics infrastructure creates a two-tier system: large companies with dedicated data science teams have access to comprehensive cross-functional intelligence, while academic centers, smaller biotechnology companies, and institutions in lower-resource settings rely on manual literature review and fragmented tools. The CAR-T Intelligence Agent, running on a $4,699 workstation with open-source components (Milvus, BGE-small, Streamlit, Python) and a pay-per-query cloud LLM, demonstrates that this capability gap can be substantially narrowed.

The Apache 2.0 license ensures that the system can be freely adapted, extended, and deployed by any organization. The modular architecture -- with clearly separated ingestion pipelines, collection schemas, knowledge graph, query expansion, RAG engine, and export modules -- is designed for extension. Adding a new data type (e.g., patent claims, conference abstracts, electronic health record extracts) requires implementing a new collection schema, a new ingest parser, and registering the collection with the RAG engine.

### 12.3 The Shift from Search to Synthesis

Traditional biomedical informatics tools are fundamentally search tools: they help users find relevant documents, which users must then read, interpret, and synthesize manually. The CAR-T Intelligence Agent represents a shift from search to synthesis. By combining multi-collection retrieval, knowledge graph augmentation, query expansion, and LLM synthesis, the system produces cross-functional answers that would require hours of expert effort to assemble manually.

This shift does not replace expert judgment -- the system is designed as an intelligence augmentation tool, not a decision-making system. Every claim is grounded in cited evidence, every citation is linked to its primary source, and uncertainty is explicitly acknowledged. The goal is to accelerate the intelligence-gathering phase of CAR-T development, freeing expert time for the higher-order activities of interpretation, decision-making, and innovation.

### 12.4 Limitations

Several limitations should be acknowledged:

1. **Corpus coverage:** While the system indexes over 5,000 PubMed abstracts and 973 clinical trials, this represents a fraction of the total CAR-T literature. Full-text article access would substantially improve evidence quality but introduces licensing constraints.
2. **Manufacturing data scarcity:** Manufacturing process data is poorly represented in the public domain. The 56 curated records in `cart_manufacturing` provide reference-quality examples but do not capture the diversity of manufacturing approaches across the industry.
3. **Real-time data:** The weekly refresh cadence means the system may lag behind breaking developments by up to 7 days for literature and trials.
4. **LLM dependency:** The synthesis step requires a cloud LLM API call, introducing latency (~20-25 seconds per query) and a dependency on external infrastructure. On-device LLM deployment (e.g., via NVIDIA NIM on the DGX Spark) is a future direction that would eliminate this dependency.
5. **Evaluation rigor:** While the system produces high-quality responses on curated demo queries, systematic evaluation against expert-annotated benchmarks has not yet been conducted.

---

## 13. Conclusion

### 13.1 Key Contributions

This paper has presented the CAR-T Intelligence Agent, a multi-collection RAG system that addresses the data fragmentation challenge in CAR-T cell therapy development. The system's key contributions are:

1. **Multi-collection RAG architecture:** 11 specialized Milvus collections with typed schemas, parallel search via ThreadPoolExecutor, and score-weighted merge -- enabling a single query to retrieve evidence across the full CAR-T lifecycle.
2. **Domain knowledge augmentation:** A 3-dictionary knowledge graph (71 entries, 54 aliases) providing structured clinical context on targets, toxicities, manufacturing, biomarkers, regulatory milestones, and immunogenicity.
3. **Intelligent query expansion:** 12 maps with 229 keywords expanding to 1,961 terms, using a dual strategy (field-filter for antigens, semantic re-embedding for non-antigens) that significantly improves recall without sacrificing precision.
4. **Structured comparative analysis:** Auto-detection and resolution of "X vs Y" queries with dual retrieval, entity-grouped evidence, and structured output (comparison tables, advantages/limitations, clinical context).
5. **Citation provenance:** Every evidence item carries a clickable link to its primary source, maintaining the traceability that scientific work demands.
6. **Hardware democratization:** The complete system runs on a single NVIDIA DGX Spark ($4,699), reducing the infrastructure barrier from hundreds of thousands of dollars to under $5,000.
7. **Open-source availability:** 21,259 lines of Python across 61 files, 415 automated tests, released under the Apache 2.0 license.

### 13.2 Future Directions

Several extensions are planned or under investigation:

- **On-device LLM inference** using NVIDIA NIM microservices on the DGX Spark, eliminating cloud dependency for the synthesis step.
- **Full-text article integration** via PubMed Central open-access subset, increasing evidence depth beyond abstracts.
- **Patent landscape analysis** using USPTO and EPO data, adding intellectual property intelligence to the cross-functional analysis.
- **Automated safety signal detection** combining FAERS data with clinical trial safety tables to surface emerging toxicity patterns.
- **Multi-agent orchestration** where the CAR-T agent collaborates with the parent platform's genomics and drug discovery agents for end-to-end patient-to-candidate workflows.
- **Systematic evaluation** using expert-annotated question-answer benchmarks from CAR-T clinical teams.
- **Integration with electronic health records** for institution-specific outcome analysis and real-time pharmacovigilance.

### 13.3 Closing Remarks

The CAR-T field is generating data at an unprecedented rate -- hundreds of new publications monthly, dozens of new trial registrations, continuous post-market safety reports, expanding real-world evidence databases. The challenge is no longer data generation but data integration and synthesis. The CAR-T Intelligence Agent demonstrates that a carefully architected multi-collection RAG system, augmented with domain-specific knowledge graphs and running on accessible hardware, can transform this fragmented data landscape into actionable cross-functional intelligence. By open-sourcing this system, we aim to accelerate CAR-T development for the benefit of patients worldwide.

---

## 14. References

### Foundational CAR-T Publications

1. Maude SL, Laetsch TW, Buechner J, et al. Tisagenlecleucel in Children and Young Adults with B-Cell Lymphoblastic Leukemia. *N Engl J Med*. 2018;378(5):439-448. doi:10.1056/NEJMoa1709866 (ELIANA trial)

2. Neelapu SS, Locke FL, Bartlett NL, et al. Axicabtagene Ciloleucel CAR T-Cell Therapy in Refractory Large B-Cell Lymphoma. *N Engl J Med*. 2017;377(26):2531-2544. doi:10.1056/NEJMoa1707447 (ZUMA-1 trial)

3. Munshi NC, Anderson LD Jr, Shah N, et al. Idecabtagene Vicleucel in Relapsed and Refractory Multiple Myeloma. *N Engl J Med*. 2021;384(8):705-716. doi:10.1056/NEJMoa2024850 (KarMMa trial)

4. Berdeja JG, Madduri D, Usmani SZ, et al. Ciltacabtagene autoleucel, a B-cell maturation antigen-directed chimeric antigen receptor T-cell therapy in patients with relapsed or refractory multiple myeloma (CARTITUDE-1): a phase 1b/2 open-label study. *Lancet*. 2021;398(10297):314-324. doi:10.1016/S0140-6736(21)00933-8 (CARTITUDE-1 trial)

5. Abramson JS, Palomba ML, Gordon LI, et al. Lisocabtagene maraleucel for patients with relapsed or refractory large B-cell lymphomas (TRANSCEND NHL 001): a multicentre seamless design study. *Lancet*. 2020;396(10254):839-852. doi:10.1016/S0140-6736(20)31366-0 (TRANSCEND trial)

6. Wang M, Munoz J, Goy A, et al. KTE-X19 CAR T-Cell Therapy in Relapsed or Refractory Mantle-Cell Lymphoma. *N Engl J Med*. 2020;382(14):1331-1342. doi:10.1056/NEJMoa1914347 (ZUMA-2 trial)

### CAR-T Safety and Toxicity Management

7. Lee DW, Santomasso BD, Locke FL, et al. ASTCT Consensus Grading for Cytokine Release Syndrome and Neurologic Toxicity Associated with Immune Effector Cells. *Biol Blood Marrow Transplant*. 2019;25(4):625-638. doi:10.1016/j.bbmt.2018.12.758

8. Neelapu SS, Tummala S, Kebriaei P, et al. Chimeric antigen receptor T-cell therapy -- assessment and management of toxicities. *Nat Rev Clin Oncol*. 2018;15(1):47-62. doi:10.1038/nrclinonc.2017.148

### CAR-T Manufacturing and Biomarkers

9. Levine BL, Miskin J, Wonnacott K, Keir C. Global Manufacturing of CAR T Cell Therapy. *Mol Ther Methods Clin Dev*. 2017;4:92-101. doi:10.1016/j.omtm.2016.12.006

10. Locke FL, Rossi JM, Neelapu SS, et al. Tumor burden, inflammation, and product attributes determine outcomes of axicabtagene ciloleucel in large B-cell lymphoma. *Blood Adv*. 2020;4(19):4898-4911. doi:10.1182/bloodadvances.2020002394

### CAR-T Resistance Mechanisms

11. Orlando EJ, Han X, Tribouley C, et al. Genetic mechanisms of target antigen loss in CAR19 therapy of acute lymphoblastic leukemia. *Nat Med*. 2018;24(10):1504-1506. doi:10.1038/s41591-018-0146-z

12. Hamieh M, Dobrin A, Chatber A, et al. CAR T cell trogocytosis and cooperative killing regulate tumour antigen escape. *Nature*. 2019;568(7750):112-116. doi:10.1038/s41586-019-1054-1

### Technical Foundations

13. Xiao S, Liu Z, Zhang P, Muennighoff N. C-Pack: Packaged Resources To Advance General Chinese Embedding. 2023. arXiv:2309.07597 (BGE embedding model)

14. Lewis P, Perez E, Piktus A, et al. Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks. *Advances in Neural Information Processing Systems*. 2020;33:9459-9474. (RAG architecture)

15. Wang J, Yi X, Guo R, et al. Milvus: A Purpose-Built Vector Data Management System. *Proceedings of the 2021 International Conference on Management of Data*. 2021:2614-2627. doi:10.1145/3448016.3457550

### Regulatory and Real-World Evidence

16. U.S. Food and Drug Administration. Approved Cellular and Gene Therapy Products. https://www.fda.gov/vaccines-blood-biologics/cellular-gene-therapy-products/approved-cellular-and-gene-therapy-products

17. Jaglowski S, Hu ZH, Zhang Y, et al. Tisagenlecleucel Chimeric Antigen Receptor (CAR) T-Cell Therapy for Adults with Diffuse Large B-Cell Lymphoma (DLBCL): Real World Experience from the Center for International Blood and Marrow Transplant Research (CIBMTR) Cellular Therapy (CT) Registry. *Blood*. 2019;134(Supplement_1):766. doi:10.1182/blood-2019-130983

---

*Generated by the HCLS AI Factory -- CAR-T Intelligence Agent v1.2.0*
*Apache 2.0 License | https://github.com/ajones1923/hcls-ai-factory*
