<!--
===========================================================================
  PROJECT BIBLE — Precision Oncology Intelligence Agent
  Single-source-of-truth technical specification
===========================================================================

  Version : 1.0.0
  Author  : Adam Jones
  Date    : March 2026
  License : Apache 2.0
  Platform: NVIDIA DGX Spark ($3,999)
===========================================================================
-->

# PROJECT BIBLE — Precision Oncology Intelligence Agent

| Field       | Value                                              |
|-------------|----------------------------------------------------|
| Version     | 1.0.0                                              |
| Author      | Adam Jones                                         |
| Date        | March 2026                                         |
| License     | Apache 2.0                                         |
| Platform    | NVIDIA DGX Spark ($3,999)                          |
| LLM         | Claude Sonnet 4.6 (Anthropic)                      |
| Embeddings  | BAAI/bge-small-en-v1.5 (384-dim)                   |
| Vector DB   | Milvus 2.4 Standalone                              |
| UI          | Streamlit (port 8526)                              |
| API         | FastAPI (port 8527)                                |

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
10. [Case Management](#10-case-management)
11. [Therapy Ranking](#11-therapy-ranking)
12. [Trial Matching](#12-trial-matching)
13. [Export System](#13-export-system)
14. [API Reference](#14-api-reference)
15. [UI Guide](#15-ui-guide)
16. [Metrics & Monitoring](#16-metrics--monitoring)
17. [Configuration](#17-configuration)
18. [Docker Deployment](#18-docker-deployment)
19. [Testing](#19-testing)
20. [Scripts & Data Seeding](#20-scripts--data-seeding)
21. [File Inventory](#21-file-inventory)
22. [Dependencies](#22-dependencies)
23. [Future Roadmap](#23-future-roadmap)

---

## 1. Executive Summary

The Precision Oncology Intelligence Agent is a RAG-powered clinical
decision-support system purpose-built for molecular tumor board (MTB)
workflows. It operates as Stage 2.5 of the HCLS AI Factory pipeline,
sitting between the genomics pipeline (Stage 1) and the drug discovery
pipeline (Stage 3) to provide real-time variant interpretation, therapy
ranking, clinical trial matching, and evidence synthesis.

### Key Metrics

| Metric                  | Value                                        |
|-------------------------|----------------------------------------------|
| Python files            | 66                                           |
| Total lines of code     | ~20,490                                      |
| Milvus collections      | 11 (10 owned + 1 read-only shared)           |
| Live vectors            | 609 (526 seed + 83 knowledge graph)          |
| Actionable gene targets | 40+                                          |
| Therapy mappings        | 30+                                          |
| Resistance mechanisms   | 12+                                          |
| Pathways                | 10+                                          |
| Biomarker panels        | 20+                                          |
| Test files / cases      | 10 files, 556 test cases (all passing)       |
| Docker services         | 6                                            |
| Seed data files         | 10 JSON files, ~773 KB total                 |
| Dependencies            | 23 packages                                  |

---

## 2. Vision & Mission

**Vision.** Every oncology patient receives an evidence-grounded,
guideline-concordant treatment recommendation within minutes of genomic
profiling -- not weeks.

**Mission.** Deliver a self-contained, GPU-accelerated precision oncology
intelligence platform that:

- Interprets somatic and germline variants with CIViC/OncoKB evidence
- Ranks therapies by evidence level with resistance awareness
- Matches patients to active clinical trials using hybrid deterministic
  and semantic search
- Generates MTB-ready packets in Markdown, JSON, PDF, and FHIR R4
- Runs entirely on a single NVIDIA DGX Spark ($3,999) with no cloud
  dependencies

---

## 3. System Overview

### Architecture Diagram (logical)

```
+--------------------------------------------------------------------+
|                     HCLS AI Factory Pipeline                       |
|                                                                    |
|  Stage 1              Stage 2.5               Stage 3              |
|  Genomics Pipeline    Oncology Intelligence   Drug Discovery       |
|  (Parabricks/         Agent                   (BioNeMo/MolMIM/    |
|   DeepVariant)        +-----------------+     DiffDock)            |
|                       | Streamlit UI    |                          |
|  FASTQ -> VCF ------->| (port 8526)     |-----> Lead Compounds    |
|                       | FastAPI         |                          |
|  genomic_evidence --->| (port 8527)     |                          |
|  (read-only)          +----|--------|---+                          |
|                            |        |                              |
|                       +----v--------v----+                         |
|                       | Milvus 2.4       |                         |
|                       | 11 collections   |                         |
|                       | 609 vectors      |                         |
|                       +------------------+                         |
+--------------------------------------------------------------------+
```

### Component Stack

| Layer          | Technology                                       |
|----------------|--------------------------------------------------|
| LLM            | Claude Sonnet 4.6 (Anthropic API)                |
| Embeddings     | BAAI/bge-small-en-v1.5, 384-dim, COSINE          |
| Vector DB      | Milvus 2.4 Standalone, IVF_FLAT index            |
| API framework  | FastAPI 0.109+, Uvicorn                          |
| UI framework   | Streamlit 1.30+                                  |
| Configuration  | Pydantic v2 BaseSettings, ONCO_ env prefix       |
| VCF parsing    | cyvcf2                                           |
| PDF export     | ReportLab 4.0+                                   |
| FHIR export    | fhir.resources 7.0+                               |
| Monitoring     | Prometheus client, custom /metrics endpoint       |
| Scheduling     | APScheduler 3.10+                                |
| Compute        | NVIDIA DGX Spark, CUDA 12.x                     |

---

## 4. Collections Catalog

All 11 collections use BGE-small-en-v1.5 embeddings (dim=384) with
IVF_FLAT indexing and COSINE similarity. Index params: nlist=1024,
search nprobe=16.

### 4.1 onco_variants (weight 0.18)

Actionable somatic/germline variants from CIViC and OncoKB.

| Field                  | Type          | Max Length |
|------------------------|---------------|------------|
| id (PK)                | VARCHAR       | 100        |
| embedding              | FLOAT_VECTOR  | 384-dim    |
| gene                   | VARCHAR       | 50         |
| variant_name           | VARCHAR       | 100        |
| variant_type           | VARCHAR       | 30         |
| cancer_type            | VARCHAR       | 50         |
| evidence_level         | VARCHAR       | 20         |
| drugs                  | VARCHAR       | 500        |
| civic_id               | VARCHAR       | 20         |
| vrs_id                 | VARCHAR       | 100        |
| text_summary           | VARCHAR       | 3000       |
| clinical_significance  | VARCHAR       | 200        |
| allele_frequency       | FLOAT         | --         |

**Seed data:** 90 records (variant_seed_data.json)

### 4.2 onco_literature (weight 0.16)

PubMed/PMC/preprint literature chunks tagged by cancer type.

| Field        | Type          | Max Length |
|--------------|---------------|------------|
| id (PK)      | VARCHAR       | 100        |
| embedding    | FLOAT_VECTOR  | 384-dim    |
| title        | VARCHAR       | 500        |
| text_chunk   | VARCHAR       | 3000       |
| source_type  | VARCHAR       | 20         |
| year         | INT64         | --         |
| cancer_type  | VARCHAR       | 50         |
| gene         | VARCHAR       | 50         |
| variant      | VARCHAR       | 100        |
| keywords     | VARCHAR       | 1000       |
| journal      | VARCHAR       | 200        |

**Seed data:** 60 records (literature_seed_data.json)

### 4.3 onco_therapies (weight 0.14)

Approved and investigational therapies with mechanism of action.

| Field                  | Type          | Max Length |
|------------------------|---------------|------------|
| id (PK)                | VARCHAR       | 100        |
| embedding              | FLOAT_VECTOR  | 384-dim    |
| drug_name              | VARCHAR       | 200        |
| category               | VARCHAR       | 30         |
| targets                | VARCHAR       | 200        |
| approved_indications   | VARCHAR       | 500        |
| resistance_mechanisms  | VARCHAR       | 500        |
| evidence_level         | VARCHAR       | 20         |
| text_summary           | VARCHAR       | 3000       |
| mechanism_of_action    | VARCHAR       | 500        |

**Seed data:** 64 records (therapy_seed_data.json)

### 4.4 onco_guidelines (weight 0.12)

NCCN/ASCO/ESMO guideline recommendations.

| Field                 | Type          | Max Length |
|-----------------------|---------------|------------|
| id (PK)               | VARCHAR       | 100        |
| embedding             | FLOAT_VECTOR  | 384-dim    |
| org                   | VARCHAR       | 20         |
| cancer_type           | VARCHAR       | 50         |
| version               | VARCHAR       | 50         |
| year                  | INT64         | --         |
| key_recommendations   | VARCHAR       | 3000       |
| text_summary          | VARCHAR       | 3000       |
| evidence_level        | VARCHAR       | 20         |

**Seed data:** 45 records (guideline_seed_data.json)

### 4.5 onco_trials (weight 0.10)

ClinicalTrials.gov summaries with biomarker eligibility criteria.

| Field               | Type          | Max Length |
|---------------------|---------------|------------|
| id (PK)             | VARCHAR       | 20         |
| embedding           | FLOAT_VECTOR  | 384-dim    |
| title               | VARCHAR       | 500        |
| text_summary        | VARCHAR       | 3000       |
| phase               | VARCHAR       | 30         |
| status              | VARCHAR       | 30         |
| sponsor             | VARCHAR       | 200        |
| cancer_types        | VARCHAR       | 200        |
| biomarker_criteria  | VARCHAR       | 500        |
| enrollment          | INT64         | --         |
| start_year          | INT64         | --         |
| outcome_summary     | VARCHAR       | 2000       |

**Seed data:** 55 records (trial_seed_data.json)

### 4.6 onco_biomarkers (weight 0.08)

Predictive and prognostic biomarkers with assay details.

| Field            | Type          | Max Length |
|------------------|---------------|------------|
| id (PK)          | VARCHAR       | 100        |
| embedding        | FLOAT_VECTOR  | 384-dim    |
| name             | VARCHAR       | 100        |
| biomarker_type   | VARCHAR       | 30         |
| cancer_types     | VARCHAR       | 200        |
| predictive_value | VARCHAR       | 500        |
| testing_method   | VARCHAR       | 500        |
| clinical_cutoff  | VARCHAR       | 500        |
| text_summary     | VARCHAR       | 3000       |
| evidence_level   | VARCHAR       | 20         |

**Seed data:** 50 records (biomarker_seed_data.json)

### 4.7 onco_resistance (weight 0.07)

Resistance mechanisms and bypass strategies for targeted therapies.

| Field                  | Type          | Max Length |
|------------------------|---------------|------------|
| id (PK)                | VARCHAR       | 100        |
| embedding              | FLOAT_VECTOR  | 384-dim    |
| primary_therapy        | VARCHAR       | 200        |
| gene                   | VARCHAR       | 50         |
| mechanism              | VARCHAR       | 500        |
| bypass_pathway         | VARCHAR       | 200        |
| alternative_therapies  | VARCHAR       | 500        |
| text_summary           | VARCHAR       | 3000       |

**Seed data:** 50 records (resistance_seed_data.json)

### 4.8 onco_pathways (weight 0.06)

Signaling pathways, cross-talk, and druggable nodes.

| Field                | Type          | Max Length |
|----------------------|---------------|------------|
| id (PK)              | VARCHAR       | 100        |
| embedding            | FLOAT_VECTOR  | 384-dim    |
| name                 | VARCHAR       | 100        |
| key_genes            | VARCHAR       | 500        |
| therapeutic_targets  | VARCHAR       | 300        |
| cross_talk           | VARCHAR       | 500        |
| text_summary         | VARCHAR       | 3000       |

**Seed data:** 35 records (pathway_seed_data.json)

### 4.9 onco_outcomes (weight 0.04)

Real-world treatment outcome records for precision oncology.

| Field                    | Type          | Max Length |
|--------------------------|---------------|------------|
| id (PK)                  | VARCHAR       | 100        |
| embedding                | FLOAT_VECTOR  | 384-dim    |
| case_id                  | VARCHAR       | 100        |
| therapy                  | VARCHAR       | 200        |
| cancer_type              | VARCHAR       | 50         |
| response                 | VARCHAR       | 20         |
| duration_months          | FLOAT         | --         |
| toxicities               | VARCHAR       | 500        |
| biomarkers_at_baseline   | VARCHAR       | 500        |
| text_summary             | VARCHAR       | 3000       |

**Seed data:** 40 records (outcome_seed_data.json)

### 4.10 onco_cases (weight 0.02)

De-identified patient case snapshots for similarity search.

| Field            | Type          | Max Length |
|------------------|---------------|------------|
| id (PK)          | VARCHAR       | 100        |
| embedding        | FLOAT_VECTOR  | 384-dim    |
| patient_id       | VARCHAR       | 100        |
| cancer_type      | VARCHAR       | 50         |
| stage            | VARCHAR       | 20         |
| variants         | VARCHAR       | 1000       |
| biomarkers       | VARCHAR       | 1000       |
| prior_therapies  | VARCHAR       | 500        |
| text_summary     | VARCHAR       | 3000       |

**Seed data:** 37 records (cases_seed_data.json)

### 4.11 genomic_evidence (weight 0.03) -- READ-ONLY

VCF-derived genomic evidence populated by Stage 1 (genomics pipeline).
The oncology agent reads from this collection but never writes to it.

| Field                  | Type          | Max Length |
|------------------------|---------------|------------|
| id (PK)                | VARCHAR       | 200        |
| embedding              | FLOAT_VECTOR  | 384-dim    |
| chrom                  | VARCHAR       | 10         |
| pos                    | INT64         | --         |
| ref                    | VARCHAR       | 500        |
| alt                    | VARCHAR       | 500        |
| qual                   | FLOAT         | --         |
| gene                   | VARCHAR       | 50         |
| consequence            | VARCHAR       | 100        |
| impact                 | VARCHAR       | 20         |
| genotype               | VARCHAR       | 10         |
| text_summary           | VARCHAR       | 2000       |
| clinical_significance  | VARCHAR       | 200        |
| rsid                   | VARCHAR       | 20         |
| disease_associations   | VARCHAR       | 500        |
| am_pathogenicity       | FLOAT         | --         |
| am_class               | VARCHAR       | 30         |

### Collection Weight Distribution

```
onco_variants      ████████████████████  0.18
onco_literature    ████████████████      0.16
onco_therapies     ██████████████        0.14
onco_guidelines    ████████████          0.12
onco_trials        ██████████            0.10
onco_biomarkers    ████████              0.08
onco_resistance    ███████               0.07
onco_pathways      ██████                0.06
onco_outcomes      ████                  0.04
genomic_evidence   ███                   0.03
onco_cases         ██                    0.02
                                    Sum: 1.00
```

---

## 5. Knowledge Graph

**Module:** `src/knowledge.py` (1,662 lines)

The knowledge graph provides curated, structured domain knowledge that is
injected into LLM prompts alongside retrieved evidence. It consists of
four primary dictionaries and one panel registry.

### 5.1 ACTIONABLE_TARGETS (~40 genes)

Each entry contains:
- `gene`, `full_name`, `cancer_types`, `key_variants`
- `targeted_therapies`, `combination_therapies`
- `resistance_mutations`, `pathway`, `evidence_level`
- `description` (free-text clinical summary)

**Representative genes:** BRAF, EGFR, ALK, ROS1, KRAS, HER2, NTRK,
RET, MET, FGFR, PIK3CA, IDH1, IDH2, BRCA1, BRCA2, TP53, PTEN,
CDKN2A, STK11, ESR1, ERBB2, NRAS, APC, VHL, KIT, PDGFRA, FLT3,
NPM1, DNMT3A, and others.

### 5.2 THERAPY_MAP (~30 drugs)

Maps drug names (lowercase) to structured metadata:
- `brand_name`, `category`, `drug_class`
- `guideline` (reference to NCCN/ESMO recommendation)

### 5.3 RESISTANCE_MAP (~12 mechanisms)

Maps drug names to documented resistance mechanisms:
- `mutation` (e.g., EGFR T790M for erlotinib)
- `next_line` (recommended subsequent therapies)
- `mechanism` (description of resistance biology)

### 5.4 PATHWAY_MAP (~10 pathways)

Maps oncogenic signaling pathways to:
- `key_genes`, `therapeutic_targets`, `cross_talk`
- Pathway descriptions for LLM context injection

### 5.5 BIOMARKER_PANELS (~20 panels)

Maps biomarker identifiers to clinical decision rules:
- `marker` name, `threshold`, `positive_values`
- `drugs` (recommended therapies), `evidence_level`
- `guideline` (reference text)

### 5.6 Helper Functions

- `lookup_gene(query)` -- Return knowledge context for gene mentions
- `lookup_therapy(query)` -- Return knowledge context for therapy mentions
- `lookup_resistance(query)` -- Return resistance mechanism context
- `lookup_pathway(query)` -- Return pathway context
- `lookup_biomarker(query)` -- Return biomarker context
- `get_target_context(gene)` -- Return full ACTIONABLE_TARGETS entry
- `classify_variant_actionability(gene, variant)` -- Return evidence tier

---

## 6. Query Expansion

**Module:** `src/query_expansion.py` (812 lines)

Domain-aware query expansion maps oncology keywords to lists of
semantically related terms to improve RAG retrieval recall. The module
defines expansion dictionaries across 12 categories:

| Category            | Example Input    | Example Expansions                                |
|---------------------|------------------|---------------------------------------------------|
| Cancer types        | NSCLC            | lung adenocarcinoma, EGFR-mutant lung, ALK-pos    |
| Genes               | EGFR             | L858R, exon 19 deletion, T790M, C797S             |
| Therapies           | osimertinib      | Tagrisso, 3rd-gen EGFR TKI, CNS penetration       |
| Biomarkers          | TMB              | tumor mutational burden, mut/Mb, KEYNOTE-158       |
| Pathways            | MAPK             | RAS-RAF-MEK-ERK, RTK signaling, BRAF cascade      |
| Resistance          | T790M            | gatekeeper mutation, osimertinib, 3rd-gen TKI      |
| Clinical terms      | PFS              | progression-free survival, median PFS, HR          |
| Trial terms         | Phase 3          | randomized, pivotal, registration trial            |
| Immunotherapy       | checkpoint       | PD-1, PD-L1, CTLA-4, anti-PD-1, pembrolizumab     |
| Surgery/radiation   | lobectomy        | surgical resection, VATS, thoracotomy              |
| Toxicity            | pneumonitis      | ILD, interstitial lung disease, immune-related     |
| Genomics            | ctDNA            | circulating tumor DNA, liquid biopsy, MRD          |

### Expansion Function

```python
def expand_query(query: str) -> List[str]:
    """Return a list of expansion terms for the given query string."""
```

The function tokenizes the input, matches against all expansion maps,
and returns the union of all matching expansions. These terms are
appended to the original query before embedding for an expanded-search
pass through the RAG engine.

---

## 7. Data Models

**Module:** `src/models.py` (538 lines)

All models use Pydantic v2 `BaseModel`. Each domain model includes a
`to_embedding_text()` method that serializes the model into a single
string suitable for BGE-small-en-v1.5 embedding.

### 7.1 Enumerations (13)

| Enum             | Values | Description                             |
|------------------|--------|-----------------------------------------|
| CancerType       | 26     | NSCLC, SCLC, BREAST, COLORECTAL, MELANOMA, PANCREATIC, OVARIAN, PROSTATE, RENAL, BLADDER, HEAD_NECK, HEPATOCELLULAR, GASTRIC, GLIOBLASTOMA, AML, CML, ALL, CLL, DLBCL, MULTIPLE_MYELOMA, CHOLANGIOCARCINOMA, ENDOMETRIAL, THYROID, MESOTHELIOMA, SARCOMA, OTHER |
| VariantType      | 7      | SNV, INDEL, CNV_AMP, CNV_DEL, FUSION, REARRANGEMENT, SV |
| EvidenceLevel    | 5      | A (FDA-approved), B (Clinical), C (Case reports), D (Preclinical), E (Computational) |
| TherapyCategory  | 9      | TARGETED, IMMUNOTHERAPY, CHEMOTHERAPY, HORMONAL, COMBINATION, RADIOTHERAPY, CELL_THERAPY, ADC, BISPECIFIC |
| TrialPhase       | 8      | Early Phase 1, Phase 1, Phase 1/2, Phase 2, Phase 2/3, Phase 3, Phase 4, N/A |
| TrialStatus      | 9      | Not yet recruiting, Recruiting, Enrolling by invitation, Active not recruiting, Suspended, Terminated, Completed, Withdrawn, Unknown status |
| ResponseCategory | 5      | CR (complete response), PR (partial response), SD (stable disease), PD (progressive disease), NE (not evaluable) |
| BiomarkerType    | 8      | PREDICTIVE, PROGNOSTIC, DIAGNOSTIC, MONITORING, RESISTANCE, PHARMACODYNAMIC, SCREENING, THERAPEUTIC_SELECTION |
| PathwayName      | 13     | MAPK, PI3K_AKT_MTOR, DDR, CELL_CYCLE, APOPTOSIS, WNT, NOTCH, HEDGEHOG, JAK_STAT, ANGIOGENESIS, HIPPO, NF_KB, TGF_BETA |
| GuidelineOrg     | 8      | NCCN, ESMO, ASCO, WHO, CAP_AMP, FDA, EMA, AACR |
| SourceType       | 4      | PUBMED, PMC, PREPRINT, MANUAL |
| (TrialPhase)     | (included above) | |
| (TrialStatus)    | (included above) | |

### 7.2 Domain Models (10)

| Model                | Key Fields                                                    |
|----------------------|---------------------------------------------------------------|
| OncologyLiterature   | id, title, text_chunk, source_type, year, cancer_type, gene, variant, keywords, journal |
| OncologyTrial        | id (NCT), title, text_summary, phase, status, sponsor, cancer_types, biomarker_criteria, enrollment, start_year, outcome_summary |
| OncologyVariant      | id, gene, variant_name, variant_type, cancer_type, evidence_level, drugs, civic_id, vrs_id, text_summary, clinical_significance, allele_frequency |
| OncologyBiomarker    | id, name, biomarker_type, cancer_types, predictive_value, testing_method, clinical_cutoff, text_summary, evidence_level |
| OncologyTherapy      | id, drug_name, category, targets, approved_indications, resistance_mechanisms, evidence_level, text_summary, mechanism_of_action |
| OncologyPathway      | id, name (PathwayName), key_genes, therapeutic_targets, cross_talk, text_summary |
| OncologyGuideline    | id, org (GuidelineOrg), cancer_type, version, year, key_recommendations, text_summary, evidence_level |
| ResistanceMechanism  | id, primary_therapy, gene, mechanism, bypass_pathway, alternative_therapies, text_summary |
| OutcomeRecord        | id, case_id, therapy, cancer_type, response (ResponseCategory), duration_months, toxicities, biomarkers_at_baseline, text_summary |
| CaseSnapshot         | case_id (alias: id), patient_id, cancer_type, stage, variants, biomarkers, prior_therapies, text_summary, created_at, updated_at |

### 7.3 Search & Agent Models (4)

| Model                 | Purpose                                          |
|-----------------------|--------------------------------------------------|
| SearchHit             | Single Milvus search result with collection, id, score, text, metadata, label, citation, relevance |
| CrossCollectionResult | Aggregated hits across collections with hit_count, hits_by_collection() |
| ComparativeResult     | Side-by-side evidence for two entities (entity_a, entity_b) |
| MTBPacket             | Molecular Tumor Board packet: variant_table, evidence_table, therapy_ranking, trial_matches, open_questions, citations |

### 7.4 Agent I/O Models (2)

| Model         | Fields                                                     |
|---------------|------------------------------------------------------------|
| AgentQuery    | question, cancer_type?, gene?, include_genomic (bool)      |
| AgentResponse | question, answer, evidence (CrossCollectionResult), knowledge_used, timestamp, plan, report |

---

## 8. RAG Engine

**Module:** `src/rag_engine.py` (908 lines)
**Class:** `OncoRAGEngine`

### 8.1 Architecture

```
User Question
      |
      v
+---[Query Embedding]--- BGE-small-en-v1.5 with instruction prefix
      |
      v
+---[Parallel Collection Search]--- ThreadPoolExecutor (max 8 workers)
      |                              Each collection searched with per-collection
      |                              weight, filter_field, year_field
      |
      v
+---[Query Expansion Search]--- (optional) expand_query() -> re-embed -> search
      |
      v
+---[Merge & Rank]--- De-duplicate by ID, weight-adjusted score, top 30
      |
      v
+---[Knowledge Injection]--- gene, therapy, resistance, pathway, biomarker
      |
      v
+---[Prompt Assembly]--- Domain Knowledge + Evidence + Question + Instructions
      |
      v
+---[LLM Synthesis]--- Claude Sonnet 4.6 -> answer with citations
      |
      v
AgentResponse
```

### 8.2 Collection Weights

Weights are applied multiplicatively to raw cosine similarity scores:

```python
COLLECTION_CONFIG = {
    "onco_variants":    {"weight": 0.18, "label": "Variant",    "filter_field": "gene"},
    "onco_literature":  {"weight": 0.16, "label": "Literature", "filter_field": "gene",  "year_field": "year"},
    "onco_therapies":   {"weight": 0.14, "label": "Therapy"},
    "onco_guidelines":  {"weight": 0.12, "label": "Guideline",  "year_field": "year"},
    "onco_trials":      {"weight": 0.10, "label": "Trial",      "year_field": "start_year"},
    "onco_biomarkers":  {"weight": 0.08, "label": "Biomarker"},
    "onco_resistance":  {"weight": 0.07, "label": "Resistance", "filter_field": "gene"},
    "onco_pathways":    {"weight": 0.06, "label": "Pathway"},
    "onco_outcomes":    {"weight": 0.04, "label": "Outcome"},
    "onco_cases":       {"weight": 0.02, "label": "Case"},
    "genomic_evidence": {"weight": 0.03, "label": "Genomic"},
}
```

### 8.3 System Prompt

The system prompt defines the agent persona as an **Oncology Intelligence
Agent** with core competencies in:

1. Molecular profiling (TMB, MSI, CNV, fusions)
2. Variant interpretation (CIViC/OncoKB evidence levels, AMP/ASCO/CAP)
3. Therapy selection (NCCN/ESMO guideline-concordant)
4. Clinical trial matching (ClinicalTrials.gov, basket/umbrella)
5. Resistance mechanisms (on-target, bypass, lineage plasticity)
6. Biomarker assessment (TMB, MSI, PD-L1, HRD, companion diagnostics)
7. Outcomes monitoring (RECIST, survival, MRD, ctDNA)
8. Cross-modal integration (genomic-imaging-drug discovery)

Behavioral instructions enforce citation, cross-functional thinking,
resistance flagging, guideline referencing, and uncertainty acknowledgment.

### 8.4 Comparative Retrieval

The engine detects comparative questions via regex
(`compare|vs|versus|difference between|head.to.head`) and routes them to
a dual-entity retrieval pipeline:

1. Parse entity A and entity B from the question
2. Retrieve evidence independently for each entity
3. Identify shared/head-to-head evidence (intersection by ID)
4. Build a structured comparison prompt with 8 comparison axes
5. Generate comparative synthesis via LLM

### 8.5 Citation Formatting

- PubMed IDs: `[PubMed 12345](https://pubmed.ncbi.nlm.nih.gov/12345/)`
- NCT IDs: `[NCT01234567](https://clinicaltrials.gov/study/NCT01234567)`
- All others: `[Label: record_id]`

### 8.6 Relevance Classification

| Score Range | Classification |
|-------------|----------------|
| >= 0.85     | high           |
| >= 0.65     | medium         |
| < 0.65      | low            |

---

## 9. Agent Architecture

**Module:** `src/agent.py` (553 lines)
**Class:** `OncoIntelligenceAgent`

### 9.1 Plan-Search-Evaluate-Synthesize Loop

```
Question
    |
    v
[1. PLAN] -----> SearchPlan
    |               - identified_topics
    |               - target_genes (from KNOWN_GENES, 30 genes)
    |               - relevant_cancer_types (from KNOWN_CANCER_TYPES + aliases)
    |               - search_strategy: broad | targeted | comparative
    |               - sub_questions (decomposed from complex queries)
    |
    v
[2. SEARCH] ----> Cross-collection retrieval
    |               - Primary query + all sub-questions
    |               - Via OncoRAGEngine.cross_collection_search()
    |
    v
[3. EVALUATE] --> "sufficient" | "partial" | "insufficient"
    |               - sufficient: >= 3 hits from >= 2 collections
    |               - If insufficient and retries remain: broaden and retry
    |               - MAX_RETRIES = 2
    |
    v
[4. SYNTHESIZE] -> AgentResponse
                    - LLM generates answer with citations
                    - Markdown report attached
```

### 9.2 Planning Heuristics

**Gene extraction:** Scans question (uppercased) against KNOWN_GENES set
(30 genes: BRAF, EGFR, ALK, ROS1, KRAS, HER2, NTRK, RET, MET, FGFR,
PIK3CA, IDH1, IDH2, BRCA, BRCA1, BRCA2, TP53, PTEN, CDKN2A, STK11,
ESR1, ERBB2, NRAS, APC, VHL, KIT, PDGFRA, FLT3, NPM1, DNMT3A).

**Cancer type extraction:** Matches against KNOWN_CANCER_TYPES (25
canonical types) and _CANCER_ALIASES (70+ aliases mapping colloquial
names to canonical types, e.g., "lung" -> "NSCLC", "tnbc" -> "BREAST").

**Topic identification:** Keyword matching across 20+ topic triggers
(resistance, biomarker, survival, immunotherapy, clinical trial,
combination, mutation, fusion, liquid biopsy, checkpoint, TMB, MSI).

**Strategy selection:**
- `comparative` -- if question contains "compare", "vs", "versus",
  "difference between", or "head to head"
- `targeted` -- if both gene(s) and cancer type(s) are identified
- `broad` -- fallback for all other queries

**Sub-question decomposition:**
- Multiple genes -> per-gene queries
- Multiple cancer types -> per-cancer queries
- Topic-specific sub-questions (resistance, trials, biomarkers, combos)

### 9.3 Evidence Evaluation

```python
MIN_SUFFICIENT_HITS = 3
MIN_COLLECTIONS_FOR_SUFFICIENT = 2
MIN_SIMILARITY_SCORE = 0.30
```

Evidence items with scores below 0.30 are discarded. The verdict is:
- `sufficient` -- >= 3 quality hits from >= 2 distinct collections
- `partial` -- some hits but below threshold
- `insufficient` -- zero quality hits

### 9.4 Fallback Queries

When evidence is insufficient, the agent generates broader fallback
queries:
- `"{gene} oncology therapeutic implications"`
- `"{gene} mutation clinical significance"`
- `"{cancer_type} precision medicine current landscape"`

---

## 10. Case Management

**Module:** `src/case_manager.py` (516 lines)
**Class:** `OncologyCaseManager`

### 10.1 Case Creation Flow

```
Patient Data (ID, cancer type, stage, VCF/variants, biomarkers, prior therapies)
    |
    v
[VCF Parsing] -- If input is raw VCF text, parse via cyvcf2
    |              Extract: gene, variant, chrom, pos, ref, alt, consequence
    |
    v
[Variant Annotation] -- Cross-reference ACTIONABLE_TARGETS
    |                    Classify actionability (A/B/C/D/E/VUS)
    |
    v
[CaseSnapshot Creation] -- Generate UUID, build text_summary
    |                       Embed via BGE-small-en-v1.5
    |
    v
[Persist to onco_cases] -- Insert into Milvus collection
    |
    v
CaseSnapshot returned
```

### 10.2 MTB Packet Generation

The MTB (Molecular Tumor Board) packet is a structured decision-support
document containing:

| Section           | Content                                          |
|-------------------|--------------------------------------------------|
| patient_summary   | Cancer type, stage, variant summary              |
| variant_table     | Gene, variant, type, evidence level, drugs       |
| evidence_table    | Retrieved evidence from RAG engine               |
| therapy_ranking   | Ranked therapies (via TherapyRanker)             |
| trial_matches     | Matched clinical trials (via TrialMatcher)       |
| open_questions    | Gaps in evidence, recommended further testing    |
| citations         | Formatted PubMed/NCT citations                  |

### 10.3 VCF Parsing

**Module:** `src/utils/vcf_parser.py`

Parses VCF files using cyvcf2, extracting:
- Chromosome, position, ref/alt alleles
- Quality score, genotype
- Gene symbol (from INFO/ANN field)
- Consequence (from INFO/ANN field)
- Allele frequency (from INFO/AF or FORMAT/AD)

---

## 11. Therapy Ranking

**Module:** `src/therapy_ranker.py` (748 lines)
**Class:** `TherapyRanker`

### 11.1 Ranking Algorithm

```
Patient Profile (cancer_type, variants, biomarkers, prior_therapies)
    |
    v
[Step 1: Variant-Driven Therapies]
    |   ACTIONABLE_TARGETS lookup for each gene/variant
    |   Evidence level from knowledge graph
    |
    v
[Step 2: Biomarker-Driven Therapies]
    |   MSI-H -> pembrolizumab, nivolumab, dostarlimab (Level A)
    |   TMB-H (>=10 mut/Mb) -> pembrolizumab (Level A)
    |   HRD/BRCA -> olaparib, rucaparib, niraparib, talazoparib (Level A/B)
    |   PD-L1 TPS >=50% -> pembrolizumab first-line (Level A)
    |   NTRK fusion -> larotrectinib, entrectinib (Level A)
    |   PTEN loss -> alpelisib (Level C)
    |   + BIOMARKER_PANELS registry check
    |
    v
[Step 3: Evidence Level Sort]
    |   A (FDA-approved) > B (Clinical) > C (Case reports) > D > E
    |
    v
[Step 4: Resistance Check]
    |   RESISTANCE_MAP: mutation-level resistance
    |   _DRUG_CLASS_GROUPS: same-mechanism class resistance
    |       egfr_tki_1g, egfr_tki_2g, egfr_tki_3g, alk_tki,
    |       braf_inhibitor, mek_inhibitor, anti_pd1, anti_pdl1,
    |       parp_inhibitor, kras_g12c
    |
    v
[Step 5: Contraindication Check]
    |   Same drug previously used -> flag
    |   Same drug_class as prior failed therapy -> flag
    |
    v
[Step 6: Supporting Evidence Retrieval]
    |   Search onco_therapies + onco_literature for each drug
    |
    v
[Step 6.5: Combination Therapy Identification]
    |   Known FDA-approved combos:
    |     dabrafenib + trametinib (BRAF, COMBI-d/v)
    |     encorafenib + binimetinib (BRAF, COLUMBUS)
    |     encorafenib + cetuximab (BRAF CRC, BEACON)
    |     ipilimumab + nivolumab (dual checkpoint)
    |     lenvatinib + pembrolizumab (endometrial/RCC)
    |     trastuzumab + pertuzumab (HER2, CLEOPATRA)
    |
    v
[Step 7: Final Ranking]
    Clean therapies first (sorted by evidence level)
    Flagged therapies after (resistance/contraindication)
    Assign rank 1..N
```

### 11.2 Output Schema

Each therapy dict contains:

```
rank, drug_name, brand_name, category, targets, evidence_level,
supporting_evidence, resistance_flag, resistance_detail,
contraindication_flag, guideline_recommendation, source,
source_gene/source_variant/source_biomarker
```

---

## 12. Trial Matching

**Module:** `src/trial_matcher.py` (513 lines)
**Class:** `TrialMatcher`

### 12.1 Hybrid Matching Strategy

```
Patient Profile (cancer_type, biomarkers, stage, age)
    |
    v
[Step 1: Deterministic Filter]
    |   Cancer type (fuzzy via 18+ alias groups)
    |   Open statuses: Recruiting, Active, Enrolling by invitation,
    |                  Not yet recruiting
    |   Query onco_trials with Milvus filter expressions
    |
    v
[Step 2: Semantic Search]
    |   Build eligibility query: "{cancer_type} clinical trial stage {X} {biomarkers}"
    |   Embed and search onco_trials via vector similarity
    |
    v
[Step 3: Merge & Deduplicate]
    |   Union by trial_id, keep best score from either source
    |
    v
[Step 4: Composite Scoring]
    |   biomarker_match = 0.40 * (fraction of patient biomarkers in criteria)
    |   semantic_score  = 0.25 * (vector similarity score)
    |   phase_weight    = 0.20 * (Phase 3=1.0, Phase 2/3=0.9, Phase 2=0.8, ...)
    |   status_weight   = 0.15 * (Recruiting=1.0, Active=0.6, Enrolling=0.8, ...)
    |   composite = (biomarker + semantic + phase + status) * age_penalty
    |
    v
[Step 5: Explanation Generation]
    |   matched_criteria: biomarkers found in trial criteria
    |   unmatched_criteria: biomarkers not confirmed
    |   explanation: human-readable match rationale
    |
    v
Ranked trial list (top_k, default 10)
```

### 12.2 Scoring Weights

| Component        | Weight | Range  |
|------------------|--------|--------|
| Biomarker match  | 0.40   | 0.0-1.0 |
| Semantic score   | 0.25   | 0.0-1.0 |
| Phase weight     | 0.20   | 0.5-1.0 |
| Status weight    | 0.15   | 0.3-1.0 |

**Age penalty:** 1.0 (in range or unknown) or 0.5 (out of range).
Detected via regex patterns: "Age >= 18", "18-75 years", etc.

---

## 13. Export System

**Module:** `src/export.py` (1,055 lines)

Four export formats, each accepting dict, MTBPacket, or string input.

### 13.1 Markdown Export

**Function:** `export_markdown(mtb_packet_or_response, title=None) -> str`

Sections: Header (patient/meta), Clinical Summary, Somatic Variant Profile
(table), Biomarker Summary, Evidence Summary, Therapy Ranking (table),
Clinical Trial Matches, Pathway Context, Known Resistance Mechanisms,
Open Questions, Disclaimer.

### 13.2 JSON Export

**Function:** `export_json(mtb_packet_or_response) -> dict`

Standardized schema with `meta` block (format, version, generated_at,
pipeline, author) plus all clinical sections.

### 13.3 PDF Export

**Function:** `export_pdf(mtb_packet_or_response, output_path) -> str`

NVIDIA-themed PDF via ReportLab with:
- Green header bar (RGB 118, 185, 0) with white title
- Structured tables for variants, therapies, trials
- Alternating row colors (whitesmoke/white)
- Footer disclaimer in gray 7pt text

**Convenience:** `markdown_to_pdf(markdown_text) -> bytes` for
converting any Markdown string to PDF bytes.

### 13.4 FHIR R4 Export

**Function:** `export_fhir_r4(mtb_packet, patient_id) -> dict`

Generates a FHIR R4 Bundle (type=collection) containing:

| Resource            | Content                                     |
|---------------------|---------------------------------------------|
| Patient             | Identifier with urn:hcls-ai-factory:patient |
| Observation (N)     | One per variant (LOINC 69548-6)             |
| Observation (TMB)   | Tumor mutation burden (LOINC 94076-7)       |
| Observation (MSI)   | Microsatellite instability (LOINC 81695-9)  |
| Specimen            | Tumor tissue (SNOMED 119376003)             |
| Condition           | Cancer diagnosis (SNOMED-coded)             |
| MedicationRequest   | Therapy recommendations (top 10)            |
| DiagnosticReport    | Master genomic report (LOINC 81247-9)       |

**SNOMED codes** for 22 cancer types are provided in
`FHIR_SNOMED_CANCER_CODES`.

---

## 14. API Reference

**Module:** `api/main.py` + `api/routes/` (5 route modules)
**Base URL:** `http://localhost:8527`
**Docs:** `http://localhost:8527/docs` (Swagger UI)

### 14.1 Core Endpoints (in main.py)

| Method | Path             | Description                              |
|--------|------------------|------------------------------------------|
| GET    | /                | Service info (docs URL, health URL)      |
| GET    | /health          | Health check with collection stats       |
| GET    | /collections     | List all collections with entity counts  |
| POST   | /query           | Full RAG query (retrieve + LLM answer)   |
| POST   | /search          | Evidence-only vector search (no LLM)     |
| POST   | /find-related    | Cross-collection entity linking          |
| GET    | /knowledge/stats | Knowledge base aggregate statistics      |
| GET    | /metrics         | Prometheus-compatible metrics             |

#### POST /query

```json
Request:
{
  "question": "What are the treatment options for EGFR L858R NSCLC?",
  "cancer_type": "nsclc",
  "gene": "EGFR",
  "top_k": 10
}

Response:
{
  "question": "...",
  "answer": "LLM-generated answer with citations...",
  "processing_time_ms": 1234.5
}
```

#### POST /search

```json
Request:
{
  "question": "BRAF V600E resistance mechanisms",
  "top_k": 10
}

Response:
{
  "results": [{"collection": "...", "id": "...", "score": 0.92, "text": "..."}],
  "count": 10,
  "processing_time_ms": 45.2
}
```

#### GET /health

```json
Response:
{
  "status": "healthy",
  "collections": {"onco_variants": 90, "onco_literature": 60, ...},
  "total_vectors": 609,
  "version": "0.1.0",
  "services": {
    "milvus": true,
    "embedder": true,
    "rag_engine": true,
    "intelligence_agent": true,
    "case_manager": true,
    "trial_matcher": true,
    "therapy_ranker": true
  }
}
```

### 14.2 Meta-Agent Router (api/routes/meta_agent.py)

| Method | Path          | Description                            |
|--------|---------------|----------------------------------------|
| POST   | /api/ask      | Unified Q&A with comparative detection |

### 14.3 Cases Router (api/routes/cases.py)

| Method | Path                  | Description                      |
|--------|-----------------------|----------------------------------|
| POST   | /api/cases            | Create patient case              |
| GET    | /api/cases/{id}       | Retrieve case by ID              |
| POST   | /api/cases/{id}/mtb   | Generate MTB packet              |
| GET    | /api/cases/{id}/variants | List case variants            |

### 14.4 Trials Router (api/routes/trials.py)

| Method | Path                        | Description                      |
|--------|-----------------------------|----------------------------------|
| POST   | /api/trials/match           | Match clinical trials            |
| POST   | /api/trials/match-case/{id} | Match trials for existing case   |
| POST   | /api/therapies/rank         | Rank therapies for patient       |

### 14.5 Reports Router (api/routes/reports.py)

| Method | Path                          | Description                |
|--------|-------------------------------|----------------------------|
| POST   | /api/reports/generate         | Generate report            |
| GET    | /api/reports/{id}/{format}    | Export (markdown/json/pdf/fhir) |

### 14.6 Events Router (api/routes/events.py)

| Method | Path          | Description                            |
|--------|---------------|----------------------------------------|
| GET    | /api/events           | Event log with pagination              |
| GET    | /api/events/{event_id}| Retrieve specific event                |

### 14.7 CORS Configuration

Allowed origins (configurable via `ONCO_CORS_ORIGINS`):
- `http://localhost:8080` (Landing page)
- `http://localhost:8526` (Streamlit UI)
- `http://localhost:8527` (API self)

### 14.8 Request Size Limit

Default: 10 MB (configurable via `ONCO_MAX_REQUEST_SIZE_MB`).
Returns HTTP 413 if exceeded.

---

## 15. UI Guide

**Module:** `app/oncology_ui.py`
**URL:** `http://localhost:8526`

The Streamlit UI provides five tabs for clinical decision support:

### Tab 1: Case Workbench

- **Create Case:** Enter patient ID, cancer type, stage, biomarkers,
  prior therapies. Upload VCF file or paste VCF text.
- **Manage Cases:** View existing cases, case details, variant tables.
- **Generate MTB Packet:** One-click MTB packet generation for any case.
  Includes variant table, evidence, therapy ranking, trial matches.

### Tab 2: Evidence Explorer

- **Free-text Q&A:** Enter any oncology question. The agent plans,
  searches, evaluates, and synthesizes a cited answer.
- **Comparative mode:** Questions with "vs" or "compare" trigger
  side-by-side evidence retrieval.
- **Citation links:** Clickable PubMed and ClinicalTrials.gov links.

### Tab 3: Trial Finder

- **Patient profile input:** Cancer type, stage, biomarkers, age.
- **Composite scoring:** Each trial shows match_score, matched_criteria,
  unmatched_criteria, and a natural-language explanation.
- **Phase/status filters:** Filter by trial phase or recruitment status.

### Tab 4: Therapy Ranker

- **Input:** Cancer type, variants (gene + variant), biomarkers, prior
  therapies.
- **Output:** Ranked therapy table with evidence level, resistance flags,
  contraindication flags, and supporting evidence citations.
- **Combination detection:** Automatically identifies FDA-approved combo
  regimens when component drugs are present.

### Tab 5: Outcomes Dashboard

- **Knowledge stats:** Actionable targets count, therapy map size,
  resistance mechanisms count, pathway count.
- **Collection metrics:** Vector counts per collection, total vectors.
- **Event log:** Recent API events and system activity.

---

## 16. Metrics & Monitoring

**Module:** `src/metrics.py` (362 lines)

### 16.1 Prometheus Metrics

The agent exposes Prometheus-compatible metrics at `GET /metrics` in text
format. Metrics include:

| Metric                         | Type      | Description                   |
|--------------------------------|-----------|-------------------------------|
| onco_agent_up                  | Gauge     | Service availability (0/1)    |
| onco_collection_vectors        | Gauge     | Vector count per collection   |
| onco_query_duration_seconds    | Histogram | RAG query latency             |
| onco_search_duration_seconds   | Histogram | Vector search latency         |
| onco_embedding_duration_seconds| Histogram | Embedding generation latency  |
| onco_llm_tokens_total          | Counter   | LLM token consumption         |
| onco_queries_total             | Counter   | Total queries processed       |
| onco_errors_total              | Counter   | Error count by type           |

### 16.2 No-Op Fallback

When `prometheus_client` is not installed, the metrics module provides
silent no-op stubs (`_NoOpMetric`) so the rest of the codebase can call
metric methods without guarding imports.

### 16.3 Health Endpoint

`GET /health` returns structured JSON with:
- `status`: "healthy" or "degraded"
- `collections`: per-collection vector counts
- `total_vectors`: aggregate count
- `version`: API version
- `services`: boolean status for each component (milvus, embedder,
  rag_engine, intelligence_agent, case_manager, trial_matcher,
  therapy_ranker)

---

## 17. Configuration

**Module:** `config/settings.py` (134 lines)
**Class:** `OncoSettings(BaseSettings)` -- Pydantic v2

All settings are overridable via environment variables with the `ONCO_`
prefix. The `.env` file is loaded automatically.

### 17.1 Complete Settings Reference

| Setting                        | Type   | Default                            | Description                              |
|--------------------------------|--------|------------------------------------|------------------------------------------|
| PROJECT_ROOT                   | Path   | (auto-detected)                    | Project root directory                   |
| DATA_DIR                       | Path   | PROJECT_ROOT/data                  | Data directory                           |
| CACHE_DIR                      | Path   | PROJECT_ROOT/cache                 | Cache directory                          |
| REFERENCE_DIR                  | Path   | PROJECT_ROOT/reference             | Reference data directory                 |
| RAG_PIPELINE_ROOT              | Path   | (auto-detected)                    | Path to rag-chat-pipeline                |
| MILVUS_HOST                    | str    | localhost                          | Milvus server hostname                   |
| MILVUS_PORT                    | int    | 19530                              | Milvus gRPC port                         |
| COLLECTION_LITERATURE          | str    | onco_literature                    | Literature collection name               |
| COLLECTION_TRIALS              | str    | onco_trials                        | Trials collection name                   |
| COLLECTION_VARIANTS            | str    | onco_variants                      | Variants collection name                 |
| COLLECTION_BIOMARKERS          | str    | onco_biomarkers                    | Biomarkers collection name               |
| COLLECTION_THERAPIES           | str    | onco_therapies                     | Therapies collection name                |
| COLLECTION_PATHWAYS            | str    | onco_pathways                      | Pathways collection name                 |
| COLLECTION_GUIDELINES          | str    | onco_guidelines                    | Guidelines collection name               |
| COLLECTION_RESISTANCE          | str    | onco_resistance                    | Resistance collection name               |
| COLLECTION_OUTCOMES            | str    | onco_outcomes                      | Outcomes collection name                 |
| COLLECTION_CASES               | str    | onco_cases                         | Cases collection name                    |
| COLLECTION_GENOMIC             | str    | genomic_evidence                   | Read-only genomic collection             |
| EMBEDDING_MODEL                | str    | BAAI/bge-small-en-v1.5             | Sentence-transformer model name          |
| EMBEDDING_DIM                  | int    | 384                                | Embedding vector dimension               |
| EMBEDDING_BATCH_SIZE           | int    | 32                                 | Batch size for embedding generation      |
| LLM_PROVIDER                   | str    | anthropic                          | LLM provider                             |
| LLM_MODEL                      | str    | claude-sonnet-4-20250514           | LLM model identifier                     |
| ANTHROPIC_API_KEY              | str?   | None                               | Anthropic API key                        |
| TOP_K                          | int    | 5                                  | Default per-collection hit limit         |
| SCORE_THRESHOLD                | float  | 0.4                                | Minimum similarity score                 |
| WEIGHT_VARIANTS                | float  | 0.18                               | onco_variants weight                     |
| WEIGHT_LITERATURE              | float  | 0.16                               | onco_literature weight                   |
| WEIGHT_THERAPIES               | float  | 0.14                               | onco_therapies weight                    |
| WEIGHT_GUIDELINES              | float  | 0.12                               | onco_guidelines weight                   |
| WEIGHT_TRIALS                  | float  | 0.10                               | onco_trials weight                       |
| WEIGHT_BIOMARKERS              | float  | 0.08                               | onco_biomarkers weight                   |
| WEIGHT_RESISTANCE              | float  | 0.07                               | onco_resistance weight                   |
| WEIGHT_PATHWAYS                | float  | 0.06                               | onco_pathways weight                     |
| WEIGHT_OUTCOMES                | float  | 0.04                               | onco_outcomes weight                     |
| WEIGHT_CASES                   | float  | 0.02                               | onco_cases weight                        |
| WEIGHT_GENOMIC                 | float  | 0.03                               | genomic_evidence weight                  |
| NCBI_API_KEY                   | str?   | None                               | NCBI/PubMed API key                      |
| PUBMED_MAX_RESULTS             | int    | 5000                               | Max PubMed ingest results                |
| CT_GOV_BASE_URL                | str    | https://clinicaltrials.gov/api/v2  | ClinicalTrials.gov API base              |
| CIVIC_BASE_URL                 | str    | https://civicdb.org/api            | CIViC API base URL                       |
| API_HOST                       | str    | 0.0.0.0                            | FastAPI bind host                        |
| API_PORT                       | int    | 8527                               | FastAPI bind port                        |
| API_BASE_URL                   | str    | http://localhost:8527              | API base URL for internal refs           |
| STREAMLIT_PORT                 | int    | 8526                               | Streamlit bind port                      |
| METRICS_ENABLED                | bool   | True                               | Enable Prometheus metrics                |
| SCHEDULER_INTERVAL             | str    | 168h                               | Background refresh interval              |
| CONVERSATION_MEMORY_DEPTH      | int    | 3                                  | Turns of conversation context            |
| CITATION_STRONG_THRESHOLD      | float  | 0.75                               | Score threshold for strong citation      |
| CITATION_MODERATE_THRESHOLD    | float  | 0.60                               | Score threshold for moderate citation    |
| CROSS_MODAL_ENABLED            | bool   | True                               | Enable cross-modal triggers              |
| CROSS_MODAL_THRESHOLD          | float  | 0.40                               | Cross-modal trigger threshold            |
| GENOMIC_TOP_K                  | int    | 5                                  | Genomic evidence retrieval limit         |
| IMAGING_TOP_K                  | int    | 5                                  | Imaging evidence retrieval limit         |
| TRIAL_WEIGHT_BIOMARKER         | float  | 0.40                               | Trial matching biomarker weight          |
| TRIAL_WEIGHT_SEMANTIC          | float  | 0.25                               | Trial matching semantic weight           |
| TRIAL_WEIGHT_PHASE             | float  | 0.20                               | Trial matching phase weight              |
| TRIAL_WEIGHT_STATUS            | float  | 0.15                               | Trial matching status weight             |
| MIN_SUFFICIENT_HITS            | int    | 3                                  | Agent: min hits for sufficient           |
| MIN_COLLECTIONS_FOR_SUFFICIENT | int    | 2                                  | Agent: min collections for sufficient    |
| MIN_SIMILARITY_SCORE           | float  | 0.30                               | Agent: min score for quality evidence    |
| CORS_ORIGINS                   | str    | localhost:8080,8526,8527           | Comma-separated CORS origins             |
| MAX_REQUEST_SIZE_MB            | int    | 10                                 | Max HTTP request body size               |
| PDF_BRAND_COLOR_R              | int    | 118                                | PDF header bar red component             |
| PDF_BRAND_COLOR_G              | int    | 185                                | PDF header bar green component           |
| PDF_BRAND_COLOR_B              | int    | 0                                  | PDF header bar blue component            |

---

## 18. Docker Deployment

**File:** `docker-compose.yml`

### 18.1 Services (6)

| Service            | Image / Build     | Port(s)      | Role                            |
|--------------------|-------------------|--------------|---------------------------------|
| milvus-etcd        | quay.io/coreos/etcd:v3.5.5 | (internal) | Milvus metadata store       |
| milvus-minio       | minio/minio:2023-03-20 | (internal) | Milvus object storage        |
| milvus-standalone  | milvusdb/milvus:v2.4-latest | 19530, 9091 | Vector database           |
| onco-streamlit     | (local build)     | 8526         | Clinical UI                     |
| onco-api           | (local build)     | 8527         | REST API (uvicorn, 2 workers)  |
| onco-setup         | (local build)     | (none)       | One-shot: create collections + seed |

### 18.2 Volumes

| Volume       | Mount Point        | Purpose                  |
|--------------|--------------------|--------------------------|
| etcd_data    | /etcd              | etcd persistent state    |
| minio_data   | /minio_data        | MinIO object storage     |
| milvus_data  | /var/lib/milvus    | Milvus vector data       |

### 18.3 Network

All services join the `onco-network` bridge network.

### 18.4 Startup Sequence

```
1. milvus-etcd starts (healthcheck: etcdctl endpoint health)
2. milvus-minio starts (healthcheck: curl minio/health/live)
3. milvus-standalone starts after etcd + minio are healthy
   (healthcheck: curl localhost:9091/healthz, start_period: 60s)
4. onco-setup runs once: create collections + seed all 10 data files
5. onco-streamlit and onco-api start after Milvus is healthy
```

### 18.5 Quick Start

```bash
cp .env.example .env           # Add ANTHROPIC_API_KEY
docker compose up -d           # Start all 6 services
docker compose logs -f onco-setup  # Watch seed progress
# Access UI at http://localhost:8526
# Access API at http://localhost:8527/docs
```

---

## 19. Testing

**Directory:** `tests/` (10 files, 4,584 lines, 556 test cases)

All tests use pytest with fixtures defined in `conftest.py`.

### 19.1 Test File Inventory

| File                    | Lines | Focus                                    |
|-------------------------|-------|------------------------------------------|
| test_models.py          | 644   | All 13 enums, 14 Pydantic models, to_embedding_text, validation |
| test_integration.py     | 785   | 4 patient profiles (NSCLC-EGFR, CRC-KRAS, Melanoma-BRAF, Breast-HER2), full pipeline end-to-end |
| test_knowledge.py       | 603   | ACTIONABLE_TARGETS completeness, THERAPY_MAP, RESISTANCE_MAP, PATHWAY_MAP, BIOMARKER_PANELS, helper functions |
| test_export.py          | 439   | Markdown structure, JSON schema, PDF generation, FHIR R4 bundle validity |
| test_therapy_ranker.py  | 363   | Variant-driven ranking, biomarker-driven ranking, resistance detection, contraindication flags, combo identification |
| test_trial_matcher.py   | 363   | Deterministic search, semantic search, composite scoring, age penalty, explanation generation |
| test_case_manager.py    | 332   | Case creation, VCF parsing, MTB packet generation, case retrieval |
| test_rag_engine.py      | 301   | Multi-collection search, weighting, citation formatting, comparative detection, relevance classification |
| test_collections.py     | 276   | Schema validation, collection creation, insert/search operations |
| test_agent.py           | 264   | Planning (gene/cancer extraction, strategy selection), evidence evaluation, report generation |

### 19.2 Running Tests

```bash
# All tests
pytest tests/ -v

# Specific module
pytest tests/test_therapy_ranker.py -v

# With coverage
pytest tests/ --cov=src --cov-report=html
```

---

## 20. Scripts & Data Seeding

**Directory:** `scripts/` (18 files)

### 20.1 Setup Scripts

| Script                    | Purpose                                      |
|---------------------------|----------------------------------------------|
| setup_collections.py      | Create all 11 Milvus collections (--drop-existing) |

### 20.2 Seed Scripts (10)

Each seed script reads its corresponding JSON file from
`data/reference/` and inserts records into the appropriate Milvus
collection with BGE-small-en-v1.5 embeddings.

| Script               | Collection       | Records | Source File                  |
|----------------------|------------------|---------|------------------------------|
| seed_variants.py     | onco_variants    | 90      | variant_seed_data.json       |
| seed_therapies.py    | onco_therapies   | 64      | therapy_seed_data.json       |
| seed_literature.py   | onco_literature  | 60      | literature_seed_data.json    |
| seed_trials.py       | onco_trials      | 55      | trial_seed_data.json         |
| seed_biomarkers.py   | onco_biomarkers  | 50      | biomarker_seed_data.json     |
| seed_resistance.py   | onco_resistance  | 50      | resistance_seed_data.json    |
| seed_guidelines.py   | onco_guidelines  | 45      | guideline_seed_data.json     |
| seed_outcomes.py     | onco_outcomes    | 40      | outcome_seed_data.json       |
| seed_cases.py        | onco_cases       | 37      | cases_seed_data.json         |
| seed_pathways.py     | onco_pathways    | 35      | pathway_seed_data.json       |
| seed_knowledge.py    | (knowledge graph)| 83      | (programmatic)               |

**Total seed records:** 526 (collections) + 83 (knowledge graph) = 609

### 20.3 Ingest Scripts (3)

| Script                      | Source                    | Purpose                  |
|-----------------------------|---------------------------|--------------------------|
| ingest_civic.py             | CIViC API                 | Actionable variant ingest |
| ingest_pubmed.py            | PubMed/NCBI E-Utilities   | Literature ingest        |
| ingest_clinical_trials.py   | ClinicalTrials.gov v2 API | Trial ingest             |

### 20.4 Validation Scripts (2)

| Script                | Purpose                                        |
|-----------------------|------------------------------------------------|
| test_rag_pipeline.py  | End-to-end RAG pipeline smoke test             |
| validate_e2e.py       | Full system validation with test queries       |

### 20.5 Ingest Parsers (src/ingest/)

| Parser                     | Lines | Source Format                    |
|----------------------------|-------|----------------------------------|
| base.py                    | --    | Abstract base class              |
| civic_parser.py            | --    | CIViC JSON API response          |
| oncokb_parser.py           | --    | OncoKB API response              |
| literature_parser.py       | --    | PubMed XML                       |
| clinical_trials_parser.py  | --    | ClinicalTrials.gov JSON          |
| guideline_parser.py        | --    | Guideline structured data        |
| pathway_parser.py          | --    | Pathway structured data          |
| resistance_parser.py       | --    | Resistance mechanism data        |
| outcome_parser.py          | --    | Treatment outcome data           |

**Total ingest module:** 9 parsers + base class, 1,793 lines.

---

## 21. File Inventory

### 21.1 Source Code (src/)

| File                     | Lines | Description                               |
|--------------------------|-------|-------------------------------------------|
| knowledge.py             | 1,662 | Knowledge graph: targets, therapies, resistance, pathways, biomarkers |
| export.py                | 1,055 | Markdown, JSON, PDF, FHIR R4 export       |
| rag_engine.py            | 908   | Multi-collection RAG engine                |
| query_expansion.py       | 812   | Domain-aware query expansion maps          |
| therapy_ranker.py        | 748   | Evidence-based therapy ranking             |
| collections.py           | 665   | 11 collection schemas + OncoCollectionManager |
| agent.py                 | 553   | Plan-search-evaluate-synthesize agent      |
| models.py                | 538   | 14 Pydantic models, 13 enums              |
| case_manager.py          | 516   | VCF parsing, case creation, MTB packets   |
| trial_matcher.py         | 513   | Hybrid trial matching                      |
| cross_modal.py           | 383   | Cross-modal trigger (genomic/imaging)      |
| metrics.py               | 362   | Prometheus metrics + no-op stubs           |
| scheduler.py             | --    | Background task scheduler                  |
| ingest/ (10 files)       | 1,793 | 9 parsers + base class                     |
| utils/ (3 files)         | 668   | VCF parser, PubMed client                  |
| **Total src/**           | **~11,440** |                                       |

### 21.2 Tests (tests/)

| File                    | Lines |
|-------------------------|-------|
| test_integration.py     | 785   |
| test_models.py          | 644   |
| test_knowledge.py       | 603   |
| test_export.py          | 439   |
| test_therapy_ranker.py  | 363   |
| test_trial_matcher.py   | 363   |
| test_case_manager.py    | 332   |
| test_rag_engine.py      | 301   |
| test_collections.py     | 276   |
| test_agent.py           | 264   |
| **Total tests/**        | **~4,584** |

### 21.3 API (api/)

| File                   | Lines | Description                        |
|------------------------|-------|------------------------------------|
| main.py                | 411   | FastAPI app, lifespan, core routes |
| routes/meta_agent.py   | --    | /api/ask unified Q&A               |
| routes/cases.py        | --    | /api/cases CRUD + MTB              |
| routes/trials.py       | --    | /api/trials/match, therapies/rank  |
| routes/reports.py      | --    | /api/reports/generate, export      |
| routes/events.py       | --    | /api/events log                    |

### 21.4 UI (app/)

| File            | Lines | Description              |
|-----------------|-------|--------------------------|
| oncology_ui.py  | --    | Streamlit 5-tab clinical UI |

### 21.5 Configuration (config/)

| File          | Lines | Description                    |
|---------------|-------|--------------------------------|
| settings.py   | 134   | Pydantic v2 BaseSettings       |

### 21.6 Top-Level Files

| File               | Description                            |
|--------------------|----------------------------------------|
| docker-compose.yml | 6-service Docker stack                 |
| Dockerfile         | Container build definition             |
| requirements.txt   | 57 Python dependencies                 |
| README.md          | Project README                         |
| LICENSE            | Apache 2.0                             |

### 21.7 Summary

| Category     | Files | Lines    |
|--------------|-------|----------|
| Source (src/) | ~25  | ~11,440  |
| Tests         | 10   | ~4,584   |
| API           | 6    | ~700     |
| UI            | 1    | ~500     |
| Config        | 1    | 134      |
| Scripts       | 18   | ~2,500   |
| Data (seed)   | 10   | ~773 KB  |
| Infra         | 3    | ~350     |
| **Total**     | **~66** | **~20,490** |

---

## 22. Dependencies

**File:** `requirements.txt` (23 packages)

### 22.1 Core Framework

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| pydantic            | >= 2.0    | Data validation, models        |
| pydantic-settings   | >= 2.7    | Environment-based configuration |
| loguru              | >= 0.7.0  | Structured logging             |

### 22.2 Vector Database

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| pymilvus            | >= 2.4.0  | Milvus Python client           |

### 22.3 Embeddings

| Package               | Version   | Purpose                      |
|------------------------|-----------|------------------------------|
| sentence-transformers  | >= 2.2.0  | BGE-small-en-v1.5 model     |

### 22.4 LLM

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| anthropic           | >= 0.18.0 | Claude API client              |

### 22.5 Web / API

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| streamlit           | >= 1.30.0 | Clinical UI framework          |
| fastapi             | >= 0.109.0| REST API framework             |
| uvicorn[standard]   | >= 0.27.0 | ASGI server                    |
| python-multipart    | >= 0.0.6  | File upload support            |

### 22.6 Data Ingest

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| requests            | >= 2.31.0 | HTTP client                    |
| lxml                | >= 5.0.0  | XML parsing (PubMed)          |
| biopython           | >= 1.83   | Biological data utilities      |

### 22.7 VCF Parsing

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| cyvcf2              | >= 0.30.0 | High-performance VCF parser    |

### 22.8 Scheduling

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| apscheduler         | >= 3.10.0 | Background task scheduler      |

### 22.9 Monitoring

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| prometheus-client   | >= 0.20.0 | Prometheus metrics             |

### 22.10 Observability

| Package              | Version   | Purpose                       |
|----------------------|-----------|-------------------------------|
| opentelemetry-api    | >= 1.29.0 | Distributed tracing API      |
| opentelemetry-sdk    | >= 1.29.0 | Distributed tracing SDK      |

### 22.11 Export / Reporting

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| reportlab           | >= 4.0.0  | PDF generation                 |
| fhir.resources      | >= 7.0.0  | FHIR R4 resource models        |

### 22.12 Utilities

| Package             | Version   | Purpose                        |
|---------------------|-----------|--------------------------------|
| numpy               | >= 1.24.0 | Numerical operations           |
| tqdm                | >= 4.65.0 | Progress bars                  |
| python-dotenv       | >= 1.0.0  | .env file loading              |

---

## 23. Future Roadmap

### Phase 1: Foundation (COMPLETE)
- [x] 11 Milvus collection schemas with IVF_FLAT / COSINE indexing
- [x] Knowledge graph: 40+ actionable targets, 30+ therapies, 12+ resistance
- [x] Multi-collection RAG engine with weighted scoring
- [x] Plan-search-evaluate-synthesize agent loop
- [x] Therapy ranking with resistance and contraindication awareness
- [x] Hybrid trial matching (deterministic + semantic)
- [x] 4-format export system (Markdown, JSON, PDF, FHIR R4)
- [x] FastAPI + Streamlit dual interface
- [x] 556 passing test cases across 10 test files
- [x] Docker Compose 6-service deployment

### Phase 2: Scale & Enrich
- [ ] Live CIViC/OncoKB ingestion pipeline (scheduled refresh)
- [ ] PubMed streaming ingest (APScheduler, configurable interval)
- [ ] ClinicalTrials.gov real-time sync
- [ ] Milvus collection partitioning by cancer type
- [ ] Embedding model upgrade evaluation (BGE-large, MedCPT)
- [ ] Multi-turn conversation memory (session-scoped context)

### Phase 3: Clinical Integration
- [ ] HL7 FHIR R4 inbound integration (receive genomic reports)
- [ ] EHR-compatible CDS Hooks interface
- [ ] Audit logging and provenance tracking (21 CFR Part 11 readiness)
- [ ] Role-based access control (oncologist vs. pathologist vs. researcher)
- [ ] PDF report co-signing workflow

### Phase 4: Advanced Intelligence
- [ ] Cross-modal pipeline triggers (genomic -> imaging -> drug discovery)
- [ ] Resistance trajectory prediction (longitudinal ctDNA tracking)
- [ ] Combination therapy optimization (synergy scoring)
- [ ] Real-world evidence integration (de-identified outcomes registry)
- [ ] Multi-agent collaboration (oncology <-> biomarker <-> CAR-T agents)

### Phase 5: Enterprise & Governance
- [ ] Multi-tenant deployment with namespace isolation
- [ ] Grafana dashboard integration (pre-built oncology panels)
- [ ] Model card and bias documentation
- [ ] External validation against AACR GENIE dataset
- [ ] SOC 2 Type II compliance documentation

---

*This document is the single source of truth for the Precision Oncology
Intelligence Agent. It is generated from verified codebase inspection
and should be updated whenever the system architecture changes.*

*Generated: March 2026 | HCLS AI Factory | Apache 2.0*
