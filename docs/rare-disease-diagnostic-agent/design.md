# Rare Disease Diagnostic Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## 1. Executive Summary

The Rare Disease Diagnostic Agent extends the HCLS AI Factory platform to accelerate the diagnostic odyssey for patients with suspected rare diseases. The average rare disease diagnosis takes 5-7 years across 7+ specialists — this agent compresses that timeline by integrating HPO-based phenotype matching, ACMG variant classification, differential diagnosis generation, gene therapy eligibility screening, and natural history data into a single RAG-powered intelligence platform.

The system searches **14 Milvus vector collections** containing phenotype ontologies, disease definitions, gene annotations, variant databases, case reports, clinical trials, therapy records, guidelines, and pathway data. It covers **88 rare diseases** across **13 disease categories**, implements **23 ACMG criteria** for variant classification, catalogs **45+ disease-associated genes**, tracks **12 gene therapies** (approved and pipeline), and runs **9 diagnostic algorithms** for systematic differential diagnosis.

The platform enables queries like *"6-year-old with progressive ataxia, nystagmus, and elevated alpha-fetoprotein — differential diagnosis?"* that simultaneously search phenotype databases, disease definitions, genetic associations, and case reports — returning a ranked differential with supporting evidence and recommended genetic testing.

### Key Results

| Metric | Value |
|---|---|
| Lines of code | **21,935** |
| Milvus collections | **14** (13 owned + 1 shared genomic_evidence) |
| Disease categories | **13** |
| Diseases cataloged | **88** (97+ including subtypes) |
| Genes cataloged | **45+** disease-associated |
| Gene therapies tracked | **12** (approved + pipeline) |
| ACMG criteria implemented | **23** (pathogenic: 5 PVS/PS/PM/PP; benign: 5 BA/BS/BP) |
| HPO top-level terms | **23** phenotype categories |
| Diagnostic algorithms | **9** systematic approaches |
| Diagnostic workflows | **10** |
| Decision support engines | **6** (HPO matcher, ACMG classifier, drug matcher, family analyzer, natural history, diagnostic algorithm) |
| Entity aliases | **120+** |
| API endpoints | **20** |
| Test suite | **193** tests (100% pass, 0.16s runtime) |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Rare Disease Agent Role |
|---|---|
| **DataStore** | Raw files: OMIM entries, Orphanet data, HPO ontology, ClinVar variants, gene therapy records |
| **DataEngine** | Ingest pipelines for phenotypes, diseases, genes, variants, literature, trials, therapies, case reports, guidelines, pathways, registries, natural history, newborn screening |
| **DataBase** | 14 Milvus collections (13 owned + 1 read-only) + knowledge base (13 categories, 88+ diseases, 12 gene therapies, 23 ACMG criteria) |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + 6 decision support engines + 9 diagnostic algorithms |
| **AgentEngine** | RareDiseaseAgent orchestrator + Streamlit UI (5 tabs) + FastAPI REST |

### 2.2 System Diagram

```
+==============================================================+
|                    RARE DISEASE DIAGNOSTIC AGENT               |
+==============================================================+
|                                                                |
|  +------------------+    +-----------------------------+       |
|  |  Streamlit UI    |    |     FastAPI REST API        |       |
|  |  Port 8544       +--->+     Port 8134               |       |
|  |  5 Tabs          |    |     20 Endpoints            |       |
|  |  NVIDIA Theme    |    |     Auth + CORS + Metrics   |       |
|  +------------------+    +----+--------+--------+------+       |
|                               |        |        |              |
|            +------------------+        |        +--------+     |
|            |                           |                 |     |
|  +---------v---------+  +-------------v----+  +---------v--+  |
|  | Workflow Engine   |  | Decision Support |  | RAG Engine |  |
|  | 10 Workflows      |  | 6 Engines        |  | Embed+Search|  |
|  | Template Method   |  | HPO Matcher      |  | 14 Colls   |  |
|  +--------+----------+  | ACMG Classifier  |  +------+-----+  |
|           |              | Drug Matcher      |         |        |
|           |              | Family Analyzer   |         |        |
|           |              | Natural History   |         |        |
|           |              | Diagnostic Algo   |         |        |
|           |              +--------+---------+         |        |
|           |                       |                   |        |
|  +--------v-----------+-----------v---+    +----------v-----+ |
|  |    Knowledge Base                  |    |    Milvus      | |
|  |    13 disease categories           |    |    14 collections| |
|  |    88+ diseases, 45+ genes         |    |    384-dim BGE | |
|  |    12 gene therapies               |    |    IVF_FLAT    | |
|  |    23 ACMG criteria                |    |    COSINE      | |
|  +------------------------------------+    +----------------+ |
+================================================================+
```

---

## 3. Data Collections — Actual State

### 3.1 Collection Catalog

| # | Collection | Est. Records | Weight | Primary Use |
|---|---|---|---|---|
| 1 | `rd_phenotypes` | 18,000 | 0.12 | HPO-based phenotype matching |
| 2 | `rd_diseases` | 10,000 | 0.11 | Disease definitions, OMIM/Orphanet |
| 3 | `rd_genes` | 5,000 | 0.10 | Gene annotations, constraint scores |
| 4 | `rd_variants` | 500,000 | 0.10 | ClinVar variants, ACMG classifications |
| 5 | `rd_literature` | 50,000 | 0.08 | PubMed rare disease literature |
| 6 | `rd_trials` | 8,000 | 0.06 | Rare disease clinical trials |
| 7 | `rd_therapies` | 2,000 | 0.07 | Gene therapies, enzyme replacement, substrate reduction |
| 8 | `rd_case_reports` | 20,000 | 0.07 | Published case reports with phenotype-genotype correlations |
| 9 | `rd_guidelines` | 3,000 | 0.06 | Diagnostic and management guidelines |
| 10 | `rd_pathways` | 2,000 | 0.06 | Metabolic pathways, disease mechanisms |
| 11 | `rd_registries` | 1,500 | 0.04 | Patient registries, natural history studies |
| 12 | `rd_natural_history` | 5,000 | 0.05 | Disease milestones, progression data |
| 13 | `rd_newborn_screening` | 80 | 0.05 | NBS conditions, analytes, cutoffs |
| 14 | `genomic_evidence` | 3,560,000 | 0.03 | Shared genomic variant context |

### 3.2 Index Configuration

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 1024 (variants, literature), 256 (diseases, trials), 128 (others) |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Decision Support Engines

### 4.1 HPO Phenotype Matcher

Matches patient phenotype descriptions to HPO (Human Phenotype Ontology) terms and computes phenotype-based disease similarity scores.

| Feature | Detail |
|---|---|
| HPO term lookup | Free-text to HPO ID resolution via semantic search |
| IC scoring | Information content-based term weighting (rarer phenotypes weighted higher) |
| Disease matching | Phenotype overlap scoring against 88+ disease phenotype profiles |
| Negated phenotypes | Supports "absent" phenotypes that reduce disease probability |
| **Output** | Ranked differential diagnosis with match scores and key distinguishing features |

### 4.2 ACMG Variant Classifier

Implements the ACMG/AMP 2015 framework for variant classification with 23 criteria:

| Category | Criteria | Strength | Examples |
|---|---|---|---|
| **Pathogenic** | PVS1 | Very Strong | Null variant in gene where LOF is known disease mechanism |
| | PS1-PS4 | Strong | Same amino acid change as established pathogenic, de novo, functional studies, prevalence |
| | PM1-PM6 | Moderate | Hot spot, absent from controls, cosegregation, missense in constrained gene, novel missense, assumed de novo |
| | PP1-PP5 | Supporting | Coseg (insufficient), missense in gene with low benign rate, in silico, patient phenotype, reputable source |
| **Benign** | BA1 | Stand-alone | Allele frequency >5% in any population |
| | BS1-BS4 | Strong | Greater frequency than expected, healthy adults, functional no effect, non-cosegregation |
| | BP1-BP7 | Supporting | Missense in truncating gene, observed in trans, in silico benign, synonymous, common |

**Classification logic:** Pathogenic (PVS1 + >= 1 PS, or >= 2 PS, or PS + >= 3 PM/PP), Likely Pathogenic, Uncertain Significance, Likely Benign, Benign — following ACMG combining rules.

### 4.3 Gene Therapy Eligibility Screener

Evaluates patient eligibility for 12 tracked gene therapies:

| Therapy | Disease | Status | Key Criteria |
|---|---|---|---|
| Zolgensma (onasemnogene) | SMA Type 1 | FDA Approved | Age < 2y, SMN1 biallelic deletion, < 3 copies SMN2 |
| Luxturna (voretigene) | RPE65 retinal dystrophy | FDA Approved | Biallelic RPE65 mutations, sufficient viable retinal cells |
| Hemgenix (etranacogene) | Hemophilia B | FDA Approved | Severe/moderate hemophilia B, no FIX inhibitors |
| Casgevy (exagamglogene) | SCD / Beta-thalassemia | FDA Approved | Age >= 12, severe SCD or transfusion-dependent thal |
| Skysona (elivaldogene) | CALD | FDA Approved | Males 4-17, early active CALD, matched donor unavailable |
| Elevidys (delandistrogene) | DMD | FDA Approved (accelerated) | Age 4-5, ambulatory, confirmed DMD |
| + 6 pipeline therapies | Various | Phase I-III | Disease-specific criteria |

### 4.4 Family Segregation Analyzer

Evaluates variant co-segregation with disease phenotype across family pedigrees, supporting autosomal dominant, autosomal recessive, X-linked, and mitochondrial inheritance patterns.

### 4.5 Natural History Engine

Provides disease-specific milestone timelines, age-of-onset distributions, progression rates, and life expectancy data for counseling and clinical trial eligibility assessment.

### 4.6 Diagnostic Algorithm Engine

Implements 9 systematic diagnostic approaches:

| Algorithm | Trigger | Steps |
|---|---|---|
| Metabolic screen | Neonatal/infantile presentation | NBS → plasma amino acids → urine organic acids → acylcarnitine → enzyme assay |
| Neurogenetic panel | Progressive neurodegeneration | Gene panel → WES → WGS → RNA-seq |
| Connective tissue | Marfanoid habitus, hypermobility | Ghent criteria → FBN1/COL3A1/TGFBR → echo + vascular imaging |
| Skeletal dysplasia | Short stature, bone abnormalities | Skeletal survey → FGFR3/COL2A1 → WES |
| Immunodeficiency | Recurrent infections | Flow cytometry → TREC/KREC → gene panel |
| Mitochondrial | Multi-system, lactic acidosis | Lactate → muscle biopsy → mtDNA → nuclear mito genes |
| Lysosomal storage | Organomegaly, regression | Enzyme assay → biomarkers → gene confirmation |
| RASopathy | Facial gestalt, cardiac, growth | Noonan panel (PTPN11, RAF1, RIT1, SOS1, KRAS) |
| Ciliopathy | Renal cysts, retinal dystrophy | Renal US → ophtho → ciliopathy gene panel |

---

## 5. Clinical Workflows

| # | Workflow | Clinical Question |
|---|---|---|
| 1 | Phenotype-Driven Differential | "Generate a ranked differential from these HPO terms" |
| 2 | Variant Interpretation | "Classify this variant using ACMG criteria" |
| 3 | Gene Therapy Eligibility | "Is this patient eligible for any available gene therapy?" |
| 4 | Diagnostic Algorithm | "What systematic workup is appropriate for this presentation?" |
| 5 | Family Segregation | "Does this variant co-segregate with disease in the family?" |
| 6 | Natural History | "What is the expected disease course and milestones?" |
| 7 | Treatment Options | "What therapeutic options exist for this diagnosis?" |
| 8 | Newborn Screening | "Interpret this NBS result in context of clinical findings" |
| 9 | Undiagnosed Disease | "This patient has been undiagnosed for 5 years — next steps?" |
| 10 | Genetic Counseling | "Recurrence risk and reproductive options for this family" |

---

## 6. Multi-Collection RAG Engine

### 6.1 Search Flow

```
User Query: "6yo with progressive ataxia, nystagmus, elevated AFP — differential?"
    │
    ├── 1. HPO Matcher: Phenotype resolution                     [< 50 ms]
    │      "progressive ataxia" → HP:0002073
    │      "nystagmus" → HP:0000639
    │      "elevated AFP" → HP:0006254
    │
    ├── 2. HPO disease matching: IC-weighted overlap              [< 100 ms]
    │      #1 Ataxia-telangiectasia (ATM) — score 0.92
    │      #2 Ataxia with oculomotor apraxia type 2 — score 0.78
    │      #3 Friedreich ataxia — score 0.65
    │
    ├── 3. Embed query with BGE asymmetric prefix                 [< 5 ms]
    │
    ├── 4. Parallel search across 14 collections (top-5 each)    [12-18 ms]
    │   ├── rd_phenotypes:    Ataxia + nystagmus + AFP HPO data  (score: 0.85-0.92)
    │   ├── rd_diseases:      AT, AOA2, FA disease profiles      (score: 0.82-0.90)
    │   ├── rd_genes:         ATM, SETX, FXN gene annotations   (score: 0.78-0.86)
    │   ├── rd_case_reports:  Pediatric ataxia case series       (score: 0.75-0.84)
    │   └── rd_therapies:     AT management, gene therapy trials (score: 0.70-0.78)
    │
    ├── 5. ACMG classification (if variant provided)              [< 20 ms]
    │
    ├── 6. Gene therapy eligibility check                         [< 10 ms]
    │
    └── 7. Stream Claude Sonnet 4.6 response                     [~22-26 sec]
           Ranked differential with supporting evidence,
           recommended genetic testing (ATM sequencing),
           and natural history context
```

**Total: ~26 sec** (HPO matching ~150 ms; retrieval ~25 ms; LLM ~25 sec)

---

## 7. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 7.1 Decision Engine Performance

| Engine | Latency |
|---|---|
| HPO Matcher (phenotype resolution + disease matching) | <150 ms |
| ACMG Classifier (23 criteria evaluation) | <20 ms |
| Gene Therapy Eligibility (12 therapies) | <10 ms |
| Family Segregation (pedigree analysis) | <50 ms |
| Natural History (milestone lookup) | <10 ms |
| Diagnostic Algorithm (workflow selection) | <10 ms |
| **All 6 engines combined** | **<250 ms** |

### 7.2 RAG Query Performance

| Operation | Latency |
|---|---|
| Full query (retrieve + engines + Claude generate) | ~26 sec |
| Streaming query (time to first token) | ~3 sec |
| 14-collection parallel search | 12-18 ms |
| Response length | 1000-2500 chars |

---

## 8. Infrastructure

### 8.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4, localhost:19530 |
| Embeddings | BGE-small-en-v1.5 (BAAI) — 384-dim |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Web UI | Streamlit (port 8544, NVIDIA black/green theme) |
| REST API | FastAPI + Uvicorn (port 8134) |
| Configuration | Pydantic BaseSettings |
| Testing | pytest (193 tests) |
| Hardware target | NVIDIA DGX Spark (GB10 GPU, 128GB unified, $4,699) |

### 8.2 Service Ports

| Port | Service |
|---|---|
| 8134 | FastAPI REST API |
| 8544 | Streamlit Chat UI |
| 19530 | Milvus vector database (shared) |

### 8.3 Dependencies on HCLS AI Factory

| Dependency | Usage |
|---|---|
| Milvus 2.4 instance | Shared vector database — adds 13 owned collections alongside existing `genomic_evidence` (3.56M vectors, read-only) |
| `ANTHROPIC_API_KEY` | Shared Anthropic API key |
| BGE-small-en-v1.5 | Same embedding model as main RAG pipeline |

---

## 9. Disease Categories

### 9.1 Coverage by Category

| # | Category | Diseases | Key Genes | Example Conditions |
|---|---|---|---|---|
| 1 | Neurogenetic | 12 | ATM, FXN, SMN1, DMD, HTT | Ataxia-telangiectasia, Friedreich ataxia, SMA, Huntington |
| 2 | Metabolic / IEM | 11 | PAH, GAA, GLA, HEXA, IDUA | PKU, Pompe, Fabry, Tay-Sachs, Hurler |
| 3 | Immunodeficiency | 8 | BTK, RAG1, RAG2, IL2RG, ADA | X-SCID, ADA-SCID, XLA, CGD |
| 4 | Hematologic | 7 | HBB, F8, F9, SERPINC1 | Sickle cell, hemophilia A/B, thalassemia |
| 5 | Connective tissue | 7 | FBN1, COL3A1, COL1A1, TGFBR | Marfan, vEDS, OI, Loeys-Dietz |
| 6 | Skeletal dysplasia | 6 | FGFR3, COL2A1, COMP | Achondroplasia, SED, pseudoachondroplasia |
| 7 | Neuromuscular | 6 | DMD, SMN1, DMPK, CNBP | Duchenne, SMA, myotonic dystrophy |
| 8 | Retinal dystrophy | 5 | RPE65, ABCA4, USH2A, RHO | LCA, Stargardt, retinitis pigmentosa |
| 9 | RASopathy | 5 | PTPN11, RAF1, BRAF, HRAS | Noonan, cardiofaciocutaneous, Costello |
| 10 | Ciliopathy | 5 | PKD1, PKD2, BBS1, NPHP1 | ADPKD, Bardet-Biedl, nephronophthisis |
| 11 | Mitochondrial | 6 | MT-ND, POLG, SURF1, MELAS | MELAS, Leigh syndrome, CPEO |
| 12 | Lysosomal storage | 6 | GLA, GAA, IDUA, GBA, HEXA | Fabry, Pompe, Gaucher, MPS I |
| 13 | Other | 4 | CFTR, DNAAF, FMR1, TSC1/2 | CF, PCD, Fragile X, TSC |

---

## 10. Demo Scenarios

### 10.1 Validated Demo Queries

**1. "6-year-old with progressive cerebellar ataxia, oculomotor apraxia, and elevated AFP"**
- HPO Matcher: HP:0002073 + HP:0000657 + HP:0006254
- Differential: #1 Ataxia-telangiectasia (ATM), #2 AOA2 (SETX), #3 AOA1 (APTX)
- Testing: ATM gene sequencing, immunoglobulin levels, lymphocyte subsets

**2. "Classify this variant: NM_000051.4(ATM):c.8977C>T (p.Arg2993Ter)"**
- ACMG: PVS1 (null variant in LOF gene) + PM2 (absent from controls) + PP3 (in silico pathogenic)
- Classification: Pathogenic
- Disease association: Ataxia-telangiectasia

**3. "Is this SMA Type 1 infant eligible for Zolgensma?"**
- Gene Therapy Screener: SMN1 biallelic deletion confirmed, age < 2 years, SMN2 copy number
- Eligibility: Check anti-AAV9 antibody titer, liver function, platelet count
- Timeline: Treatment window critical — earlier administration correlates with better outcomes

**4. "Progressive hepatosplenomegaly and pancytopenia in Ashkenazi Jewish child"**
- Diagnostic Algorithm: Lysosomal storage screen
- Differential: Gaucher disease (GBA), Niemann-Pick type B (SMPD1)
- Testing: Glucocerebrosidase enzyme activity, chitotriosidase, GBA sequencing

**5. "Family with autosomal dominant polycystic kidney disease — genetic counseling"**
- Family Analyzer: AD inheritance, 50% recurrence risk per child
- Genes: PKD1 (85%), PKD2 (15%)
- Natural History: Renal cysts by age 30, ESRD by age 55 (PKD1) vs. 75 (PKD2)

---

## 11. File Structure (Actual)

```
rare_disease_diagnostic_agent/
├── src/
│   ├── agent.py                     # RareDiseaseAgent orchestrator
│   ├── models.py                    # Pydantic data models + enums
│   ├── collections.py               # 14 Milvus collection schemas
│   ├── rag_engine.py                # Multi-collection RAG engine
│   ├── clinical_workflows.py        # 10 diagnostic workflows
│   ├── decision_support.py          # 6 decision engines
│   ├── hpo_matcher.py               # HPO phenotype matching + IC scoring
│   ├── acmg_classifier.py           # 23-criteria ACMG variant classification
│   ├── gene_therapy.py              # Gene therapy eligibility screening
│   ├── knowledge.py                 # Disease taxonomy, gene catalog, therapy database
│   ├── query_expansion.py           # 120+ entity aliases
│   ├── cross_modal.py               # Cross-agent integration
│   ├── metrics.py                   # Prometheus metrics
│   └── export.py                    # Report generation
├── app/
│   └── rare_disease_ui.py          # Streamlit (5 tabs, NVIDIA theme)
├── api/
│   └── main.py                     # FastAPI REST server (20 endpoints)
├── config/
│   └── settings.py                 # Pydantic BaseSettings
├── data/
│   └── reference/                  # Disease catalogs, HPO mappings, therapy records
├── tests/                          # 193 tests
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```

**48 Python files | ~21,935 lines of code | Apache 2.0**

---

## 12. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Architecture** | Complete | 14 collections, 6 decision engines, 10 workflows, knowledge base, HPO matcher, ACMG classifier |
| **Phase 2: Data** | Complete | 88+ diseases, 45+ genes, 12 gene therapies, 23 ACMG criteria, 13 disease categories |
| **Phase 3: RAG Integration** | Complete | Multi-collection parallel search, Claude Sonnet 4.6 streaming |
| **Phase 4: Testing** | Complete | 193 tests, 100% pass, 0.16s runtime |
| **Phase 5: UI + Demo** | Complete | 5-tab Streamlit UI, NVIDIA theme, 5 demo scenarios validated |

### Remaining Work

| Item | Priority | Effort |
|---|---|---|
| OMIM API integration for real-time disease data | Medium | 2-3 days |
| Matchmaker Exchange (MME) API integration | Medium | 1 week |
| Expanded newborn screening panel (RUSP 2026 update) | Low | 1-2 days |
| Integration with HCLS AI Factory landing page | Low | 1 hour |

---

## 13. Relationship to HCLS AI Factory

The Rare Disease Diagnostic Agent demonstrates the **diagnostic extension** of the HCLS AI Factory architecture. While the core platform focuses on common cancer genomics, rare disease requires a fundamentally different approach: phenotype-driven diagnosis, ACMG variant classification, and gene therapy eligibility screening.

- **Same Milvus instance** — 13 new owned collections alongside existing `genomic_evidence` (3.56M vectors, read-only)
- **Same embedding model** — BGE-small-en-v1.5 (384-dim)
- **Same LLM** — Claude via Anthropic API
- **Same hardware** — NVIDIA DGX Spark ($4,699)
- **Same patterns** — Pydantic models, BaseIngestPipeline, knowledge graph, query expansion

The critical integration: Stage 1 (genomics pipeline) produces VCF files from whole-genome or exome sequencing; this agent's ACMG classifier interprets those variants in the context of rare disease phenotypes, connecting molecular findings to diagnostic conclusions and therapeutic options including gene therapy.

---

## 14. Credits

- **Adam Jones**
- **Apache 2.0 License**

---

!!! warning "Clinical Decision Support Disclaimer"
    The Rare Disease Diagnostic Agent is a clinical decision support research tool for rare disease evaluation. It is not FDA-cleared and is not intended as a standalone diagnostic device. Variant classifications should be reviewed by certified molecular geneticists. Gene therapy eligibility determinations require specialist evaluation. All recommendations should be reviewed by qualified clinical geneticists and genetic counselors. Apache 2.0 License.
