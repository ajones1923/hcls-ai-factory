# Rare Disease Diagnostic Agent -- Project Bible

**Version:** 1.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [Overview](#1-overview)
2. [Architecture](#2-architecture)
3. [Collections Reference](#3-collections-reference)
4. [Workflow Reference](#4-workflow-reference)
5. [API Endpoint Reference](#5-api-endpoint-reference)
6. [Knowledge Base Reference](#6-knowledge-base-reference)
7. [Decision Support Engines](#7-decision-support-engines)
8. [Query Expansion Reference](#8-query-expansion-reference)
9. [Configuration Reference](#9-configuration-reference)
10. [Port Map](#10-port-map)
11. [Tech Stack](#11-tech-stack)
12. [Data Models](#12-data-models)
13. [Cross-Agent Integration](#13-cross-agent-integration)
14. [Ingest Pipeline Reference](#14-ingest-pipeline-reference)
15. [Test Reference](#15-test-reference)

---

## 1. Overview

The Rare Disease Diagnostic Agent is an AI-powered rare disease diagnostic decision support system that integrates RAG-based evidence retrieval across 14 Milvus vector collections, 10 diagnostic workflows, 6 decision support engines, and an autonomous reasoning pipeline. It serves clinical geneticists, genetic counselors, metabolic specialists, and undiagnosed disease programs with evidence-based differential diagnosis, variant interpretation, and therapeutic matching across 13 disease categories and 97+ rare conditions.

### Key Numbers

| Metric | Value |
|---|---|
| Python files | 48 |
| Lines of code | 21,935 |
| Milvus collections | 14 |
| Diagnostic workflows | 10 |
| Decision support engines | 6 |
| API endpoints | 20 |
| Disease categories | 13 |
| Diseases cataloged | 97+ |
| Genes cataloged | 45+ |
| Gene therapies | 12 |
| ACMG criteria | 28 |
| HPO top-level terms | 23 |
| Diagnostic algorithms | 9 |
| Entity aliases | 120+ |
| Tests | 193 (100% pass, 0.16s) |
| Knowledge version | 1.0.0 |

---

## 2. Architecture

```
User --> Streamlit UI (:8544) --> FastAPI API (:8134)
                                       |
                      +----------------+----------------+
                      |                |                |
               Workflows(10)    Decision Engines(6)  RAG Engine
                      |                                |
                 Knowledge Base                  Milvus (:19530)
               (13 categories,                   14 collections
                97+ diseases,                    384-dim BGE
                12 gene therapies)               IVF_FLAT/COSINE
```

### Tiers

- **Presentation:** Streamlit (5 tabs, NVIDIA dark theme, port 8544)
- **Application:** FastAPI (20 endpoints, CORS, auth, port 8134)
- **Data:** Milvus (14 collections, BGE-small-en-v1.5 embeddings, port 19530)

---

## 3. Collections Reference

### 3.1 Full Collection Catalog

| # | Name | Est. Records | Weight | Key Fields |
|---|---|---|---|---|
| 1 | rd_phenotypes | 18,000 | 0.12 | hpo_id, name, definition, synonyms, ic_score, frequency, is_negated |
| 2 | rd_diseases | 10,000 | 0.11 | disease_id, name, omim_id, orpha_code, inheritance, prevalence, category, clinical_features |
| 3 | rd_genes | 5,000 | 0.10 | gene_symbol, gene_name, chromosome, constraint_score, disease_associations, function_summary |
| 4 | rd_variants | 500,000 | 0.10 | variant_id, gene, hgvs, classification, population_freq, clinvar_stars, review_status |
| 5 | rd_literature | 50,000 | 0.08 | pmid, title, abstract, journal, year, disease_context |
| 6 | rd_trials | 8,000 | 0.06 | nct_id, title, condition, intervention, phase, status, eligibility |
| 7 | rd_therapies | 2,000 | 0.07 | therapy_name, indication, mechanism, status, approval_year, gene_target |
| 8 | rd_case_reports | 20,000 | 0.07 | case_id, phenotypes_hpo, diagnosis, genotype, age_onset, outcome |
| 9 | rd_guidelines | 3,000 | 0.06 | guideline_id, title, source, disease, recommendation, evidence_level |
| 10 | rd_pathways | 2,000 | 0.06 | pathway_id, name, genes, enzymes, metabolites, disease_associations |
| 11 | rd_registries | 1,500 | 0.04 | registry_name, disease, organization, enrollment, country |
| 12 | rd_natural_history | 5,000 | 0.05 | disease, milestone, age_range, frequency, source |
| 13 | rd_newborn_screening | 80 | 0.05 | condition, analyte, cutoff, confirmatory_test, act_sheet |
| 14 | genomic_evidence | 3,560,000 | 0.03 | variant_id, gene, consequence, clinical_significance, evidence_text |

### 3.2 Index Configuration

All collections share: BGE-small-en-v1.5 (384-dim), IVF_FLAT, COSINE, nlist=128.

### 3.3 Workflow-Specific Weight Maps

Each of the 10 workflows has a custom weight distribution across all 14 collections. The top-weighted collection for each workflow:

| Workflow | Top Collection | Weight |
|---|---|---|
| Phenotype-Driven | rd_phenotypes | 0.22 |
| WES/WGS Interpretation | rd_variants | 0.22 |
| Metabolic Screening | rd_pathways | 0.20 |
| Dysmorphology | rd_phenotypes | 0.25 |
| Neurogenetic | rd_genes | 0.18 |
| Cardiac Genetics | rd_genes | 0.18 |
| Connective Tissue | rd_phenotypes | 0.18 |
| Inborn Errors | rd_pathways | 0.20 |
| Gene Therapy Eligibility | rd_therapies | 0.22 |
| Undiagnosed Disease | rd_phenotypes | 0.15 |

---

## 4. Workflow Reference

### 4.1 Workflow Catalog

| # | Workflow Type | Enum Value | Description |
|---|---|---|---|
| 1 | Phenotype-Driven Diagnosis | phenotype_driven | Match HPO terms to candidate diseases via BMA similarity |
| 2 | WES/WGS Interpretation | wes_wgs_interpretation | ACMG variant classification from exome/genome data |
| 3 | Metabolic Screening | metabolic_screening | Newborn screening result and metabolic profile evaluation |
| 4 | Dysmorphology | dysmorphology | Facial and skeletal feature matching for syndromic diagnosis |
| 5 | Neurogenetic Evaluation | neurogenetic | Specialized workup for neurological genetic conditions |
| 6 | Cardiac Genetics | cardiac_genetics | Cardiomyopathy and channelopathy evaluation |
| 7 | Connective Tissue | connective_tissue | Marfan, EDS, OI diagnostic workup |
| 8 | Inborn Errors | inborn_errors | Deep metabolic investigation for IEM |
| 9 | Gene Therapy Eligibility | gene_therapy_eligibility | Match patients to gene therapies and trials |
| 10 | Undiagnosed Disease | undiagnosed_disease | Multi-modal workup for unresolved cases |

### 4.2 Additional Workflow Types (API)

| Workflow Type | Enum Value | Description |
|---|---|---|
| Variant Interpretation | variant_interpretation | Focused variant-level ACMG analysis |
| Differential Diagnosis | differential_diagnosis | Broad differential with phenotype matching |
| Newborn Screening | newborn_screening | NBS result interpretation and follow-up |
| Metabolic Workup | metabolic_workup | Metabolic laboratory evaluation |
| Carrier Screening | carrier_screening | Carrier status for recessive conditions |
| Prenatal Diagnosis | prenatal_diagnosis | Prenatal genetic testing guidance |
| Natural History | natural_history | Disease course prediction |
| Therapy Selection | therapy_selection | Therapeutic option evaluation |
| Clinical Trial Matching | clinical_trial_matching | Trial eligibility assessment |
| Genetic Counseling | genetic_counseling | Counseling talking points |
| General | general | General rare disease query |

### 4.3 Workflow Output Structure

All workflows return `WorkflowResult`:

```python
WorkflowResult(
    workflow_type=DiagnosticWorkflowType,
    findings=["Finding 1", "Finding 2"],
    recommendations=["Recommend WES", "Refer to genetics"],
    guideline_references=["ACMG 2015", "GeneReviews"],
    severity=SeverityLevel.HIGH,
    cross_agent_triggers=["[CARDIOLOGY] Cardiac phenotype"],
    confidence=0.85,
    diagnostic_result=DiagnosticResult(...)
)
```

---

## 5. API Endpoint Reference

### 5.1 Core Endpoints

| Method | Path | Request | Response |
|---|---|---|---|
| GET | /health | -- | `{status, collections, vector_counts, uptime}` |
| GET | /collections | -- | `{collections: [{name, count}]}` |
| GET | /workflows | -- | `{workflows: [str]}` |
| GET | /metrics | -- | Prometheus text format |

### 5.2 Diagnostic Endpoints (v1)

| Method | Path | Request Body | Response |
|---|---|---|---|
| POST | /v1/diagnostic/query | `{query, top_k, workflow_type}` | `{response, sources, confidence}` |
| POST | /v1/diagnostic/search | `{query, collections, top_k}` | `{results: [SearchResult]}` |
| POST | /v1/diagnostic/diagnose | `PatientQuery` | `DiagnosticResult` |
| POST | /v1/diagnostic/variants/interpret | `{variants: [VariantData]}` | `{classifications: [VariantClassification]}` |
| POST | /v1/diagnostic/phenotype/match | `{hpo_terms, top_k}` | `{candidates: [DiseaseCandidate]}` |
| POST | /v1/diagnostic/therapy/search | `{disease, genotype}` | `{therapies: [TherapyMatch]}` |
| POST | /v1/diagnostic/trial/match | `{disease, genotype, criteria}` | `{trials: [TrialMatch]}` |
| POST | /v1/diagnostic/workflow/{type} | `{inputs}` | `WorkflowResult` |

### 5.3 Reference Endpoints

| Method | Path | Response |
|---|---|---|
| GET | /v1/diagnostic/disease-categories | `{categories: {name: {description, example_diseases, key_genes}}}` |
| GET | /v1/diagnostic/gene-therapies | `{therapies: [{name, disease, gene, mechanism, status}]}` |
| GET | /v1/diagnostic/acmg-criteria | `{criteria: [{code, strength, description}]}` |
| GET | /v1/diagnostic/hpo-categories | `{terms: [{hpo_id, name, description}]}` |
| GET | /v1/diagnostic/knowledge-version | `{version, last_updated, counts}` |

### 5.4 Report & Event Endpoints

| Method | Path | Purpose |
|---|---|---|
| POST | /v1/reports/generate | Generate diagnostic report (Markdown/JSON/PDF) |
| GET | /v1/reports/formats | List supported export formats |
| GET | /v1/events/stream | SSE event stream for real-time updates |

---

## 6. Knowledge Base Reference

### 6.1 Data Sources (10)

1. OMIM (Online Mendelian Inheritance in Man)
2. Orphanet Rare Disease Database
3. GeneReviews (NCBI)
4. ClinGen / ClinVar
5. Human Phenotype Ontology (HPO)
6. ACMG/AMP Standards and Guidelines (Richards et al. 2015)
7. NIH GARD (Genetic and Rare Diseases Information Center)
8. European Reference Networks (ERNs)
9. FDA Approved Cellular and Gene Therapy Products
10. Newborn Screening ACTion (ACT) Sheets -- ACMG

### 6.2 Disease Category Detail

| Category | Count | Key Genes | Diagnostic Approach |
|---|---|---|---|
| Metabolic | 28 | PAH, GBA, GLA, GAA, IDUA, ACADM, GALT, HEXA | NBS, tandem MS, enzyme assays, molecular |
| Neurological | 23 | SMN1, DMD, MECP2, HTT, FXN, SCN1A, NF1 | EMG/NCV, brain MRI, gene panels, WES/WGS |
| Hematologic | 15 | HBB, F8, F9, VWF, RPS19, FANCA | CBC, Hb electrophoresis, coagulation, breakage |
| Connective Tissue | 10 | FBN1, COL5A1, COL3A1, COL1A1, TGFBR1 | Clinical criteria, echo, molecular testing |
| Immunologic | 13 | IL2RG, JAK3, RAG1, CYBB, STAT3, WAS, BTK | Ig levels, lymphocyte subsets, TREC, DHR |
| Cancer Predisposition | 8 | TP53, MLH1, BRCA1, BRCA2, APC, RET, VHL | Pedigree, NCCN criteria, germline panels |
| Cardiac | 6+ | MYH7, MYBPC3, KCNQ1, SCN5A, PKP2 | ECG, echo, cardiac MRI, gene panels |
| Endocrine | 6+ | CYP21A2, PTPN11, FGFR1, TSHR | Hormone panels, NBS, karyotype |
| Skeletal | 6+ | FGFR3, ALPL, COL1A1, COL2A1, RUNX2 | Skeletal survey, DXA, gene panels |
| Renal | 6+ | PKD1, PKD2, COL4A5, CTNS, GLA | Renal imaging, urinalysis, gene panels |
| Pulmonary | 3+ | CFTR, SERPINA1, DNAH5 | Sweat chloride, alpha-1 AT level, ciliary EM |
| Dermatologic | 3+ | COL7A1, KRT14, XPC | Skin biopsy, gene panels |
| Ophthalmologic | 3+ | RPE65, RHO, GUCY2D | Retinal exam, ERG, gene panels |

---

## 7. Decision Support Engines

### 7.1 Engine Catalog

| # | Engine | Class | Key Algorithm | Output |
|---|---|---|---|---|
| 1 | HPO-to-Gene Matcher | HPOToGeneMatcher | BMA + IC scoring | Ranked gene list with scores |
| 2 | ACMG Variant Classifier | ACMGVariantClassifier | 28-criteria scoring | Classification + criteria + summary |
| 3 | Orphan Drug Matcher | OrphanDrugMatcher | Disease/pathway/repurposing match | TherapyMatch list |
| 4 | Diagnostic Algorithm Recommender | DiagnosticAlgorithmRecommender | 6 phenotype cluster pathways | Ordered test list with yields |
| 5 | Family Segregation Analyzer | FamilySegregationAnalyzer | Simplified LOD scoring | LOD score + ACMG evidence level |
| 6 | Natural History Predictor | NaturalHistoryPredictor | Registry-derived milestones | Milestone list with confidence |

### 7.2 Gene-HPO Association Map (14 genes)

CFTR, FBN1, SCN1A, DMD, HTT, KCNQ1, MYH7, SMN1, MECP2, PAH, GBA1, COL5A1, OTC, FGFR3

---

## 8. Query Expansion Reference

### 8.1 Entity Alias Categories

| Category | Examples | Count |
|---|---|---|
| Metabolic abbreviations | PKU, MSUD, MCAD, MPS, GSD, CDG | 20+ |
| Neurological abbreviations | SMA, DMD, CMT, TSC, NF1, HD | 10+ |
| Immunologic abbreviations | SCID, CGD, CVID, WAS, XLA | 5+ |
| Connective tissue | EDS, OI, LDS | 3+ |
| Cardiac | HCM, DCM, LQTS, ARVC | 4+ |
| Gene names | CFTR, FBN1, SMN1, etc. | 30+ |
| HPO term aliases | "seizures", "floppy baby", "short stature" | 40+ |
| **Total** | | **120+** |

---

## 9. Configuration Reference

### 9.1 Environment Variables

| Variable | Default | Description |
|---|---|---|
| RD_MILVUS_HOST | localhost | Milvus server hostname |
| RD_MILVUS_PORT | 19530 | Milvus server port |
| RD_API_PORT | 8134 | FastAPI server port |
| RD_STREAMLIT_PORT | 8544 | Streamlit UI port |
| RD_ANTHROPIC_API_KEY | (none) | Claude API key |
| RD_EMBEDDING_MODEL | BAAI/bge-small-en-v1.5 | Embedding model |
| RD_LLM_MODEL | claude-sonnet-4-6 | LLM model identifier |
| RD_SCORE_THRESHOLD | 0.4 | Minimum similarity score |
| RD_API_KEY | (empty) | API authentication key |
| RD_CORS_ORIGINS | localhost:8080,8134,8544 | Allowed CORS origins |

---

## 10. Port Map

| Port | Service | Protocol |
|---|---|---|
| 8134 | FastAPI REST API | HTTP |
| 8544 | Streamlit UI | HTTP |
| 19530 | Milvus | gRPC |
| 9091 | etcd (Milvus metadata) | gRPC |
| 9000 | MinIO (Milvus storage) | HTTP |

---

## 11. Tech Stack

| Layer | Technology |
|---|---|
| Compute | NVIDIA DGX Spark, CUDA 12.x |
| LLM | Claude (Anthropic) |
| Embeddings | BGE-small-en-v1.5 (384-dim) |
| Vector DB | Milvus 2.x (pymilvus) |
| API | FastAPI + Uvicorn |
| UI | Streamlit |
| Data Models | Pydantic v2 + pydantic-settings |
| Logging | Loguru |
| Reports | ReportLab (PDF), python-docx (DOCX) |
| Testing | pytest |
| Containerization | Docker + Docker Compose |
| Orchestration | Nextflow DSL2 (platform-level) |

---

## 12. Data Models

### 12.1 Enums

DiagnosticWorkflowType (19), InheritancePattern (7), ACMGClassification (5), VariantType (7), DiseaseCategory (14), TherapyStatus (5), Urgency (3), SeverityLevel (5), EvidenceLevel (5)

### 12.2 Pydantic Models

PatientQuery (16 fields), HPOTerm (5), DiseaseCandidate (12), VariantClassification (12), TherapyMatch (8), DiagnosticSearchResult (4), DiagnosticResult (7), WorkflowResult (8)

### 12.3 Dataclasses

SearchPlan (9 fields), CollectionConfig (6 fields)

---

## 13. Cross-Agent Integration

| Agent | Port | URL | Trigger Condition |
|---|---|---|---|
| Genomics Pipeline | 8527 | localhost:8527 | VCF data available |
| PGx Agent | 8107 | localhost:8107 | Drug metabolism variants |
| Cardiology Agent | 8126 | localhost:8126 | Cardiac phenotype detected |
| Biomarker Agent | 8529 | localhost:8529 | Biomarker stratification needed |
| Clinical Trial Agent | 8538 | localhost:8538 | Trial eligibility matching |

---

## 14. Ingest Pipeline Reference

### 14.1 Parsers

| Parser | File | Source | Target Collection |
|---|---|---|---|
| HPO Parser | ingest/hpo_parser.py | HPO OBO/JSON-LD | rd_phenotypes |
| OMIM Parser | ingest/omim_parser.py | OMIM API | rd_diseases, rd_genes |
| Orphanet Parser | ingest/orphanet_parser.py | Orphanet XML | rd_diseases |
| Gene Therapy Parser | ingest/gene_therapy_parser.py | FDA/EMA curated | rd_therapies |
| Base Ingest | ingest/base.py | -- | Abstract base class |

### 14.2 Scripts

| Script | Purpose | Command |
|---|---|---|
| setup_collections.py | Create Milvus schemas | `python scripts/setup_collections.py` |
| seed_knowledge.py | Seed knowledge base | `python scripts/seed_knowledge.py` |
| run_ingest.py | Full ingest pipeline | `python scripts/run_ingest.py` |
| generate_docx.py | MD to DOCX conversion | `python scripts/generate_docx.py` |

---

## 15. Test Reference

### 15.1 Test Summary

193 tests, 100% pass, 0.16s, 14 test files.

### 15.2 Coverage Areas

- Agent pipeline and search plan generation
- All 20 API endpoints
- 10 clinical workflows
- 14 collection schemas and weight maps
- 6 decision support engines
- End-to-end diagnostic integration
- Knowledge base data integrity
- Pydantic model validation
- Query expansion and alias resolution
- RAG engine and context assembly
- Configuration validation
- Workflow dispatch and execution

---

*Apache 2.0 License -- HCLS AI Factory*
