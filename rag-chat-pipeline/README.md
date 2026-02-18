# RAG/Chat Pipeline

[![NVIDIA](https://img.shields.io/badge/NVIDIA-DGX%20Spark-76B900?style=flat&logo=nvidia)](https://www.nvidia.com/en-us/data-center/dgx-spark/)
[![Milvus](https://img.shields.io/badge/Milvus-Vector%20DB-00A4EF?style=flat)](https://milvus.io/)
[![Claude](https://img.shields.io/badge/Anthropic-Claude-CC785C?style=flat)](https://www.anthropic.com/)
[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat&logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](../LICENSE)

**Stage 2 of the Precision Medicine to Drug Discovery AI Factory**

> Retrieval-Augmented Generation (RAG) pipeline for querying genetic variants with natural language. Transforms annotated VCF data into therapeutic intelligence using semantic search, knowledge graphs, and AI reasoning.

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                    PRECISION MEDICINE TO DRUG DISCOVERY AI FACTORY                   │
├──────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                      │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────────────┐   │
│  │  GENOMICS   │    │  RAG/CHAT   │    │   CRYO-EM   │    │ MOLECULE GENERATION │   │
│  │  PIPELINE   │───▶│  PIPELINE   │───▶│  EVIDENCE   │───▶│     (BioNeMo)       │   │
│  │             │    │ (This Repo) │    │             │    │                     │   │
│  └─────────────┘    └─────────────┘    └─────────────┘    └─────────────────────┘   │
│    FASTQ→VCF         VCF→Target        Target→Structure    Structure→Molecules      │
│    Parabricks        Milvus+Claude     PDB/EMDB            MolMIM+DiffDock          │
│                                                                                      │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Table of Contents

- [Overview](#overview)
- [From Raw Variants to Actionable Intelligence](#from-raw-variants-to-actionable-intelligence)
- [What Pharma Companies Actually Use](#what-pharma-companies-actually-use)
- [Key Features](#key-features)
- [Architecture](#architecture)
- [Annotation Pipeline](#annotation-pipeline)
  - [ClinVar: Clinical Evidence](#clinvar-clinical-evidence)
  - [AlphaMissense: AI-Predicted Pathogenicity](#alphamissense-ai-predicted-pathogenicity)
  - [VEP: Functional Consequence Prediction](#vep-functional-consequence-prediction)
- [Vector Database: Milvus](#vector-database-milvus)
- [Knowledge Connection Layer: Clinker](#knowledge-connection-layer-clinker)
- [AI Reasoning: Claude](#ai-reasoning-claude)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Demo Queries](#demo-queries)
- [Database Statistics](#database-statistics)
- [Configuration](#configuration)
- [Directory Structure](#directory-structure)
- [Troubleshooting](#troubleshooting)
- [Related Pipelines](#related-pipelines)
- [References](#references)

---

## Overview

This pipeline is the **intelligence layer** of the Precision Medicine to Drug Discovery AI Factory. It takes the VCF file from Stage 1 (Genomics Pipeline) and transforms it into actionable therapeutic insights through:

1. **Multi-source annotation** (ClinVar, AlphaMissense, VEP)
2. **Semantic search** (Milvus vector database)
3. **Knowledge connections** (Clinker - 201 genes, 100+ diseases)
4. **AI reasoning** (Claude with RAG architecture)

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         STAGE 2: RAG/CHAT PIPELINE                           │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   VCF File (11.7M variants)                                                 │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │                    ANNOTATION LAYER                           │        │
│   │                                                               │        │
│   │  ┌─────────┐    ┌──────────────┐    ┌──────────────┐        │        │
│   │  │ ClinVar │    │ AlphaMissense│    │     VEP      │        │        │
│   │  │ (Known) │    │  (AI Pred)   │    │ (Functional) │        │        │
│   │  └────┬────┘    └──────┬───────┘    └──────┬───────┘        │        │
│   │       │                │                   │                 │        │
│   │       └────────────────┼───────────────────┘                 │        │
│   │                        ▼                                     │        │
│   │              Combined Evidence                               │        │
│   │       (35,616 ClinVar + 6,831 AlphaMissense)                │        │
│   └──────────────────────────────────────────────────────────────┘        │
│                        │                                                   │
│                        ▼                                                   │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │                    VECTOR DATABASE                            │        │
│   │                                                               │        │
│   │   Milvus: 3.56M variant embeddings (BGE-small-en-v1.5)       │        │
│   │   Hybrid search: semantic + metadata filtering               │        │
│   │                                                               │        │
│   └──────────────────────────────────────────────────────────────┘        │
│                        │                                                   │
│                        ▼                                                   │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │                    CLINKER KNOWLEDGE LAYER                    │        │
│   │                                                               │        │
│   │   201 genes → Proteins → Pathways → Diseases → Drugs         │        │
│   │   Coverage: 13 therapeutic areas (85% druggable)             │        │
│   │                                                               │        │
│   └──────────────────────────────────────────────────────────────┘        │
│                        │                                                   │
│                        ▼                                                   │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │                    CLAUDE AI REASONING                        │        │
│   │                                                               │        │
│   │   RAG architecture: Grounded responses with citations         │        │
│   │   Natural language queries → Therapeutic insights            │        │
│   │                                                               │        │
│   └──────────────────────────────────────────────────────────────┘        │
│                        │                                                   │
│                        ▼                                                   │
│              Target Hypothesis → Stage 3 (Drug Discovery)                 │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## From Raw Variants to Actionable Intelligence

Once we have the VCF from the Genomics Pipeline, the next step is **annotation**—this is where genetic differences start to gain meaning.

From the roughly **11.7 million variants** identified across all chromosomes, annotation enriches each variant with biological and clinical context by linking it to:
- The **gene** it affects
- The **type of change** it causes (missense, frameshift, etc.)
- Whether it has been **observed before** in clinical databases
- **AI predictions** of pathogenicity

This process allows us to quickly separate normal human variation from the small subset of variants that may disrupt protein function or be associated with disease.

### The Filtering Funnel

```
11.7M variants (raw from VCF)
      │
      ▼ Quality filter (QUAL > 30)
3.56M high-quality variants
      │
      ▼ ClinVar annotation
35,616 clinically annotated variants
      │
      ▼ AlphaMissense prediction
6,831 AI-predicted pathogenic variants
      │
      ▼ Clinker knowledge matching
Variants in 80 druggable genes
      │
      ▼ Natural language query
Therapeutic insights for specific diseases
```

---

## What Pharma Companies Actually Use

Enterprise genomics pipelines draw from multiple tiers of annotation sources:

| Tier | Source Type | Examples | In This Pipeline |
|------|-------------|----------|------------------|
| **Clinical** | Curated databases | ClinVar, OMIM, HGMD, ClinGen | ClinVar |
| **Population** | Frequency databases | gnomAD, UK Biobank, 23andMe | (Future) |
| **AI Prediction** | Functional predictors | AlphaMissense, CADD, SpliceAI | AlphaMissense |
| **Consequence** | Functional annotation | VEP, SnpEff, ANNOVAR | VEP |

This pipeline demonstrates how these annotation layers work together:
- **ClinVar** for clinical evidence (what we know)
- **AlphaMissense** for AI-predicted pathogenicity (what AI predicts)
- **VEP** for functional consequence prediction (what the variant does)

---

## Key Features

### Multi-Source Annotation
- **ClinVar**: 35,616 clinically-annotated variants with pathogenicity classifications
- **AlphaMissense**: 6,831 AI-predicted pathogenic variants (from 71M predictions)
- **VEP**: Functional consequence annotation (missense, frameshift, splice, etc.)

### Semantic Search at Scale
- **Milvus Vector Database**: Millisecond search across 3.56M variant embeddings
- **Hybrid Search**: Combine semantic similarity with metadata filtering
- **BGE Embeddings**: State-of-the-art text embeddings (384 dimensions)

### Knowledge Graph (Clinker)
- **80 high-value genes** across 6 therapeutic areas
- **100+ disease conditions** with therapeutic connections
- **66 druggable targets** (82%) with FDA-approved drugs
- Visual knowledge paths: Variant → Gene → Protein → Pathway → Disease → Drug

### AI-Powered Reasoning
- **Claude (Anthropic)**: Advanced reasoning with RAG architecture
- **Grounded Responses**: All answers cite specific variant evidence
- **Streaming**: Real-time response generation

### File Manager
- **VCF Upload**: Upload VCF and VCF.gz files directly from the browser
- **Directory Browser**: Browse input/ and output/ directories
- **File Operations**: Download, delete, and manage files
- **File Metadata**: View size, modification date, and file type

### Therapeutic Coverage

| Therapeutic Area | Genes | Example Conditions |
|------------------|-------|-------------------|
| **Oncology** | 25 | Breast cancer, Lung cancer, Leukemia, Melanoma |
| **Neurology** | 14 | FTD, ALS, Alzheimer's, Parkinson's, Huntington's |
| **Rare Disease** | 12 | Cystic fibrosis, SMA, Muscular dystrophy, Hemophilia |
| **Cardiovascular** | 10 | Cardiomyopathy, Long QT, Hypercholesterolemia |
| **Immunology** | 8 | Rheumatoid arthritis, Psoriasis, Crohn's disease |
| **Pharmacogenomics** | 10 | Drug metabolism, Warfarin sensitivity, Chemotherapy toxicity |

---

## Architecture

### System Architecture

```
┌────────────────────────────────────────────────────────────────────────────────┐
│                              RAG/CHAT PIPELINE                                  │
├────────────────────────────────────────────────────────────────────────────────┤
│                                                                                │
│   ┌──────────────────────────────────────────────────────────────────────┐    │
│   │                        STREAMLIT UI (Port 8501)                       │    │
│   │                                                                       │    │
│   │   ┌─────────────┐  ┌──────────────────┐  ┌─────────────────────┐    │    │
│   │   │ Chat Input  │  │ Evidence Display │  │ Clinker Knowledge   │    │    │
│   │   │             │  │ (with citations) │  │ (visual graph)      │    │    │
│   │   └─────────────┘  └──────────────────┘  └─────────────────────┘    │    │
│   │                                                                       │    │
│   └───────────────────────────────────┬───────────────────────────────────┘    │
│                                       │                                        │
│   ┌───────────────────────────────────▼───────────────────────────────────┐    │
│   │                           RAG ENGINE                                   │    │
│   │                                                                       │    │
│   │   Query Analysis → Gene Expansion → Vector Search → Context Assembly │    │
│   │                                                                       │    │
│   └───────────┬───────────────────────┬───────────────────────┬───────────┘    │
│               │                       │                       │                │
│               ▼                       ▼                       ▼                │
│   ┌───────────────────┐   ┌───────────────────┐   ┌───────────────────┐       │
│   │      MILVUS       │   │     CLINKER       │   │      CLAUDE       │       │
│   │   Vector Store    │   │  Knowledge Base   │   │    LLM Client     │       │
│   │                   │   │                   │   │                   │       │
│   │  3.56M embeddings  │   │   201 genes       │   │  Anthropic API    │       │
│   │  COSINE similarity│   │   100+ diseases   │   │  Streaming SSE    │       │
│   │  IVF_FLAT index   │   │  171 drug targets │   │  RAG grounding    │       │
│   │                   │   │                   │   │                   │       │
│   └───────────────────┘   └───────────────────┘   └───────────────────┘       │
│                                                                                │
│   ┌──────────────────────────────────────────────────────────────────────┐    │
│   │                        DATA LAYER                                     │    │
│   │                                                                       │    │
│   │   ┌───────────┐   ┌──────────────┐   ┌──────────────┐               │    │
│   │   │  ClinVar  │   │ AlphaMissense│   │  VCF Parser  │               │    │
│   │   │ 4.1M vars │   │  71M preds   │   │  (cyvcf2)    │               │    │
│   │   └───────────┘   └──────────────┘   └──────────────┘               │    │
│   │                                                                       │    │
│   └──────────────────────────────────────────────────────────────────────┘    │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              DATA FLOW                                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   User Query: "What BRCA variants do I have?"                              │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 1. QUERY ANALYSIS                                             │        │
│   │    • Extract entities: BRCA → BRCA1, BRCA2                   │        │
│   │    • Identify intent: variant discovery                       │        │
│   │    • Expand genes: add related oncology genes                 │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 2. SEMANTIC SEARCH (Milvus)                                   │        │
│   │    • Embed query using BGE-small-en-v1.5                     │        │
│   │    • Search 3.56M variants by cosine similarity               │        │
│   │    • Apply metadata filter: gene IN (BRCA1, BRCA2)           │        │
│   │    • Return top-k results with scores                        │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 3. KNOWLEDGE CONNECTION (Clinker)                             │        │
│   │    • Match genes to knowledge base: BRCA1, BRCA2 → found     │        │
│   │    • Retrieve: protein, pathway, diseases, drugs             │        │
│   │    • BRCA1 → PARP inhibitors (Olaparib, Rucaparib)          │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 4. CONTEXT ASSEMBLY                                           │        │
│   │    • Format evidence with citations                          │        │
│   │    • Include Clinker knowledge connections                   │        │
│   │    • Add AlphaMissense scores where available                │        │
│   │    • Build structured prompt for Claude                      │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 5. CLAUDE REASONING                                           │        │
│   │    • Stream response with SSE                                │        │
│   │    • Explain findings in clinical context                    │        │
│   │    • Cite specific variants as evidence                      │        │
│   │    • Suggest therapeutic implications                        │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   Response: "I found 3 BRCA1/2 variants in your genome..."                │
│             + Evidence panel + Clinker visualization                       │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Annotation Pipeline

### ClinVar: Clinical Evidence

**What it is:** ClinVar is the NIH's public database of clinically interpreted genetic variants. When a clinical lab or research group determines that a variant causes disease—or confirms it's benign—they submit that interpretation to ClinVar.

**What it provides:**
- Peer-reviewed, evidence-backed classifications
- Clinical significance: Pathogenic, Likely Pathogenic, Benign, Likely Benign, VUS
- Associated disease phenotypes
- Review status indicating evidence strength
- Links to supporting publications

**Implementation:**

```python
# ClinVarAnnotator loads 4.1M GRCh38 variants at initialization
class ClinVarAnnotator:
    def __init__(self, clinvar_file: Path):
        self._variant_db = {}  # Indexed by chr_pos_ref_alt

    def annotate(self, variant: VariantEvidence) -> VariantEvidence:
        key = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
        if key in self._variant_db:
            variant.clinical_significance = data['clinical_significance']
            variant.disease_associations = data['disease_associations']
            variant.rsid = data['rsid']
        return variant
```

**Statistics:**
- Database size: 4.1 million GRCh38 variants
- Matches in HG002: 35,616 variants (1% of high-quality variants)

---

### AlphaMissense: AI-Predicted Pathogenicity

**What it is:** AlphaMissense is Google DeepMind's machine learning model that predicts whether a genetic variant will damage protein function. Built on top of AlphaFold's protein structure predictions, it asks: "Given how this protein folds in 3D, will swapping this amino acid break it?"

**Why it matters:** This allows us to assess the millions of variants that aren't in ClinVar—potential novel drug targets that haven't been clinically studied yet. AlphaFold tells us the protein's shape; AlphaMissense tells us if a mutation will break that shape.

**Implementation:**

```python
# AlphaMissenseAnnotator loads 71M predictions
class AlphaMissenseAnnotator:
    def __init__(self, alphamissense_file: Path):
        self._variant_db = {}  # 71M missense predictions

    def annotate(self, variant: VariantEvidence) -> VariantEvidence:
        key = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
        if key in self._variant_db:
            variant.am_pathogenicity = data['am_pathogenicity']  # 0.0 - 1.0
            variant.am_class = data['am_class']  # likely_benign/ambiguous/likely_pathogenic
        return variant
```

**Classification Thresholds:**
- **Likely Pathogenic**: Score > 0.564
- **Ambiguous**: Score 0.340 - 0.564
- **Likely Benign**: Score < 0.340

**Statistics:**
- Database size: 71,697,560 possible human missense variants
- Matches in HG002: 6,831 variants with pathogenicity predictions

**Novel Target Discovery:**
The combination of ClinVar (what we know) and AlphaMissense (what AI predicts) enables queries like: *"Show me high-confidence damaging variants in druggable genes that haven't been clinically studied"*—precisely the novel target discovery workflow that pharmaceutical companies use.

---

### VEP: Functional Consequence Prediction

**What it is:** VEP (Variant Effect Predictor) is Ensembl's tool for determining what type of change a variant causes. It answers: "Does this variant sit in a gene? Does it change an amino acid? Does it disrupt splicing?"

**What it provides:**
- Affected gene and transcript
- Consequence type (missense_variant, stop_gained, frameshift_variant, splice_donor_variant)
- Protein position and amino acid change
- Impact severity (HIGH, MODERATE, LOW, MODIFIER)

**How it complements other annotations:**
- **VEP** describes what the variant does structurally
- **ClinVar** provides clinical evidence of its effect
- **AlphaMissense** offers AI prediction of its impact

Together, these three annotation sources enable both clinical interpretation and novel target discovery.

---

## Vector Database: Milvus

**What it is:** Milvus is an open-source vector database purpose-built for AI applications. Traditional databases search by exact matches—"find all variants where gene equals VCP." Vector databases search by meaning.

**Why it matters:** A query about "dementia" automatically finds variants annotated with "frontotemporal lobar degeneration" or "cognitive decline" because these concepts are nearby in vector space. Researchers can ask natural questions without knowing exact terminology.

**Implementation:**

```python
class MilvusClient:
    def __init__(self, embedding_dim: int = 384):
        self.collection_name = "genomic_evidence"

    def search(self, query_embedding: np.ndarray, top_k: int = 10,
               filter_expr: Optional[str] = None) -> List[Dict]:
        # Hybrid search: semantic + metadata filtering
        results = collection.search(
            data=[query_embedding.tolist()],
            anns_field="embedding",
            param={"metric_type": "COSINE", "params": {"nprobe": 16}},
            limit=top_k,
            expr=filter_expr,  # e.g., "gene == 'BRCA1'"
            output_fields=["chrom", "pos", "gene", "clinical_significance", ...]
        )
        return results
```

**Technical Details:**
- **Embedding Model**: BGE-small-en-v1.5 (384 dimensions)
- **Index Type**: IVF_FLAT (nlist=1024)
- **Metric**: Cosine similarity
- **Collection Size**: 3.56 million variant embeddings
- **Search Latency**: <100ms

---

## Knowledge Connection Layer: Clinker

**What it is:** Clinker is the semantic layer that transforms isolated variant annotations into connected biological narratives. Annotation tells you *what* a variant is; Clinker tells you *why it matters*.

**The Connection Chain:**

```
Variant → Gene → Protein → Pathway → Disease → Drug
   │        │        │         │         │        │
   │        │        │         │         │        └── Therapeutic options
   │        │        │         │         └── Clinical relevance
   │        │        │         └── Biological context
   │        │        └── Molecular function
   │        └── Gene symbol
   └── Genomic coordinates
```

**Example Connection (VCP):**

```
rs188935092 (chr7:117559590 G>A)
      │
      ▼
Gene: VCP
      │
      ▼
Protein: p97/VCP ATPase
      │
      ▼
Pathway: Ubiquitin-proteasome system
      │
      ▼
Diseases: Frontotemporal Dementia, ALS, IBMPFD
      │
      ▼
Drugs: CB-5083 (Phase I), NMS-873, DBeQ
```

**Implementation:**

```python
KNOWLEDGE_CONNECTIONS = {
    'VCP': {
        'protein': 'p97/VCP ATPase',
        'function': 'Protein quality control, ERAD, autophagy',
        'pathway': 'Ubiquitin-proteasome system',
        'diseases': ['Frontotemporal Dementia (FTD)', 'ALS', 'IBMPFD'],
        'drugs': ['CB-5083 (Phase I)', 'NMS-873', 'DBeQ'],
        'drug_status': 'Clinical development',
        'pdb_ids': ['5FTK', '7K56', '8OOI'],
        'druggable': True,
    },
    # ... 200 more genes
}
```

**Coverage Statistics:**

| Metric | Value |
|--------|-------|
| Total Genes | 201 |
| Druggable Targets | 171 (85%) |
| Disease Conditions | 150+ |
| FDA-Approved Drugs | 100+ |
| Therapeutic Areas | 13 |

**Therapeutic Area Breakdown:**

| Area | Genes | Key Examples |
|------|-------|--------------|
| **Oncology** | 23 | BRCA1, BRCA2, EGFR, KRAS, ALK, BRAF, HER2, PD-1, PD-L1 |
| **Neurology** | 36 | VCP, GRN, C9orf72, MAPT, APP, PSEN1, LRRK2, SNCA, HTT, PINK1, TREM2, CGRP |
| **Rare Disease** | 16 | CFTR, SMN1, DMD, HBB, F8, GAA, GBA |
| **Cardiovascular** | 12 | LDLR, PCSK9, TTR, MYBPC3, SCN5A, KCNH2 |
| **Immunology** | 8 | IL6, TNF, JAK1, JAK2, IL17A, IL23A |
| **Pharmacogenomics** | 6 | CYP2D6, CYP2C19, CYP3A4, DPYD, TPMT |
| **Metabolic/Endocrine** | 22 | GLP1R, SGLT2, DPP4, PPARG, GCK, INS |
| **Infectious Disease** | 21 | HIV RT/PR/IN, HCV NS3/NS5, SARS-CoV-2 targets |
| **Respiratory** | 13 | ADRB2, IL4R, IL5, BMPR2, CFTR |
| **Ophthalmology** | 11 | VEGFA, CFH, RPE65, RHO |
| **Dermatology** | 9 | IL31RA, TYK2, IL13, COL7A1 |
| **Hematology** | 12 | SYK, THPO, F10, BTK |
| **GI/Hepatology** | 12 | ATP4A, S1PR1, THR_BETA, FXR |

---

## AI Reasoning: Claude

**What it is:** Claude is Anthropic's large language model that serves as the reasoning and communication layer. While Milvus finds relevant evidence and Clinker connects it to biological context, Claude synthesizes everything into coherent, actionable answers.

**Why it matters:** Claude doesn't hallucinate genomic facts because it's grounded in retrieved data. It acts as an expert interpreter that can explain complex genetic findings to both technical and non-technical audiences.

**RAG Architecture:**

```python
class RAGEngine:
    def query(self, user_query: str) -> Generator[str, None, None]:
        # 1. Analyze query and expand genes
        genes = self._extract_genes(user_query)
        expanded_genes = self._expand_pharmacogenomics(user_query, genes)

        # 2. Semantic search in Milvus
        query_embedding = self.embedder.embed(user_query)
        evidence = self.milvus.search(query_embedding, top_k=20)

        # 3. Get knowledge connections
        knowledge = get_knowledge_for_evidence(evidence)

        # 4. Build prompt with context
        prompt = self._build_prompt(user_query, evidence, knowledge)

        # 5. Stream Claude response
        for chunk in self.llm.stream(prompt):
            yield chunk
```

**Implementation Details:**
- **Model**: claude-sonnet-4-20250514
- **Temperature**: 0.3 (factual consistency)
- **Streaming**: Server-Sent Events (SSE)
- **Grounding**: All responses cite specific variant evidence

---

## Quick Start

### Prerequisites

- **Docker & Docker Compose** (for Milvus)
- **Python 3.10+**
- **Anthropic API Key** (for Claude)
- **VCF file** from Genomics Pipeline
- **Annotation data**: ClinVar + AlphaMissense (~2 GB)

### Download Annotation Data

```bash
# From HCLS AI Factory root (recommended):
./setup-data.sh --stage2

# This downloads ClinVar (~480 MB) and AlphaMissense (~614 MB)
# with automatic verification. See docs/DATA_SETUP.md for details.
```

### Installation

```bash
# Clone and setup
cd ~/projects/hcls-ai-factory/rag-chat-pipeline
./run.sh setup

# Configure environment
cp .env.example .env
nano .env  # Add your ANTHROPIC_API_KEY
```

### Start Services

```bash
# 1. Start Milvus vector database
./run.sh start

# 2. Ingest variants (first time only)
./run.sh ingest --annotated-only  # Fast: ~35K ClinVar variants
# OR
./run.sh ingest                   # Full: ~3.56M high-quality variants

# 3. Start chat interface
./run.sh chat
```

Open **http://localhost:8501** in your browser.

---

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd rag-chat-pipeline
```

### Step 2: Setup Virtual Environment

```bash
./run.sh setup
# OR manually:
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Step 3: Configure Environment

```bash
cp .env.example .env
```

Edit `.env`:

```bash
# Required
ANTHROPIC_API_KEY=your_key_here

# Milvus
MILVUS_HOST=localhost
MILVUS_PORT=19530

# Data paths
VCF_PATH=data/input/HG002.genome.vcf.gz
CLINVAR_PATH=data/annotations/variant_summary.txt.gz
ALPHAMISSENSE_PATH=data/annotations/AlphaMissense_hg38.tsv.gz
```

### Step 4: Start Milvus

```bash
docker-compose up -d
```

### Step 5: Download Annotation Databases

```bash
# ClinVar (automatic in ingestion script)
# AlphaMissense (614MB)
wget -P data/annotations/ https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
```

### Step 6: Ingest Variants

```bash
# Option 1: Annotated variants only (fast, demo-ready)
source venv/bin/activate && source .env
python scripts/ingest_vcf.py --annotated-only

# Option 2: Full high-quality variants (comprehensive)
python scripts/ingest_vcf.py --limit 3500000

# Option 3: With AlphaMissense
python scripts/ingest_vcf.py --alphamissense data/annotations/AlphaMissense_hg38.tsv.gz
```

### Step 7: Start Chat UI

```bash
./run.sh chat
# OR
streamlit run app/chat_ui.py --server.port 8501
```

---

## Usage

### Command Line Interface

```bash
./run.sh <command>

Commands:
  setup       Install dependencies
  start       Start Milvus database
  stop        Stop all services
  status      Check service status
  ingest      Ingest VCF into vector DB
  chat        Start Streamlit chat interface
```

### Ingestion Options

```bash
python scripts/ingest_vcf.py [OPTIONS]

Options:
  --annotated-only         Only ingest ClinVar-annotated variants
  --limit N                Maximum variants to ingest
  --drop-existing          Drop and re-create collection
  --clinvar PATH           ClinVar file path
  --alphamissense PATH     AlphaMissense file path
  --batch-size N           Embedding batch size (default: 1000)
  --use-cache              Enable embedding cache
```

### Chat Interface

The Streamlit UI provides:
- **Natural language query input**
- **Evidence panel** with expandable variant details
- **Clinker visualization** showing Gene → Protein → Disease → Drug
- **AlphaMissense scores** color-coded by pathogenicity
- **File Manager** for browsing/uploading VCF files
- **Export to Drug Discovery Pipeline**

### File Manager

The integrated File Manager (accessible via sidebar → "Files" tab) provides:

```
┌─────────────────────────────────────────────────────────────────┐
│                        FILE MANAGER                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Directory: [📁 INPUT ▾] [📁 OUTPUT]                            │
│                                                                  │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │  UPLOAD VCF FILES                                        │   │
│  │  [Choose VCF file...] (.vcf, .vcf.gz)                   │   │
│  │  [⬆️ Upload File]                                        │   │
│  └─────────────────────────────────────────────────────────┘   │
│                                                                  │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │  FILES IN INPUT                                          │   │
│  │                                                          │   │
│  │  ▶ 🧬 HG002.genome.vcf.gz                               │   │
│  │      Size: 1.2 GB | Modified: 2025-01-13 14:30          │   │
│  │      [⬇️ Download] [🗑️ Delete]                          │   │
│  │                                                          │   │
│  │  ▶ 🧬 patient_sample.vcf.gz                             │   │
│  │      Size: 856 MB | Modified: 2025-01-12 09:15          │   │
│  │      [⬇️ Download] [🗑️ Delete]                          │   │
│  │                                                          │   │
│  └─────────────────────────────────────────────────────────┘   │
│                                                                  │
│  3 files (2.1 GB) | 2 VCF files                                 │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**Features:**
- **Upload VCF/VCF.gz**: Drag and drop or browse for VCF files
- **Browse Directories**: Switch between input/ and output/ folders
- **File Details**: View size, modification date, and file type
- **Download Files**: Download any file directly to your computer
- **Delete Files**: Remove files from the pipeline

---

## Demo Queries

### Oncology

| Query | What It Does |
|-------|--------------|
| "What BRCA variants do I have?" | Finds BRCA1/BRCA2 variants, shows PARP inhibitor connections |
| "Show me lung cancer variants" | EGFR, ALK, KRAS, ROS1 variants with targeted therapy options |
| "What pathogenic variants affect cancer genes?" | Comprehensive oncology gene panel |

### Neurology

| Query | What It Does |
|-------|--------------|
| "What variants are associated with frontotemporal dementia?" | VCP, GRN, C9orf72, MAPT variants |
| "Do I have any ALS-related variants?" | SOD1, FUS, TARDBP, C9orf72 variants |
| "Find Alzheimer's disease variants" | APP, PSEN1, PSEN2, APOE variants |

### Rare Disease

| Query | What It Does |
|-------|--------------|
| "What cystic fibrosis variants do I have?" | CFTR variants with Trikafta eligibility |
| "Show me muscular dystrophy variants" | DMD variants with exon-skipping options |
| "Find sickle cell or thalassemia variants" | HBB variants with gene therapy connections |

### Cardiovascular

| Query | What It Does |
|-------|--------------|
| "What heart disease variants do I have?" | MYBPC3, MYH7, TTR, channel genes |
| "Find cholesterol-related variants" | LDLR, PCSK9, APOB variants |
| "Show me arrhythmia variants" | SCN5A, KCNH2, KCNQ1 Long QT variants |

### Pharmacogenomics

| Query | What It Does |
|-------|--------------|
| "What drug metabolism variants do I have?" | CYP2D6, CYP2C19, CYP3A4 variants |
| "Am I sensitive to warfarin?" | CYP2C9, VKORC1 variants |
| "Check for chemotherapy toxicity risk" | DPYD, TPMT, UGT1A1 variants |

---

## Database Statistics

### Current Ingestion (HG002 Full Genome)

| Metric | Value |
|--------|-------|
| **Total variants in VCF** | 11,700,000 |
| **High-quality (QUAL>30)** | 3,561,170 |
| **ClinVar annotated** | 35,616 |
| **AlphaMissense annotated** | 6,831 |
| **In Milvus vector DB** | 3,561,170 |
| **Ingestion time** | ~4 hours |

### Annotation Database Sizes

| Database | Variants | Size |
|----------|----------|------|
| ClinVar | 4,100,000 | 1.2 GB |
| AlphaMissense | 71,697,560 | 614 MB |

---

## Configuration

### Environment Variables

```bash
# API Keys
ANTHROPIC_API_KEY=sk-ant-...

# Milvus Connection
MILVUS_HOST=localhost
MILVUS_PORT=19530
MILVUS_COLLECTION=genomic_evidence

# Data Paths
VCF_PATH=data/input/HG002.genome.vcf.gz
CLINVAR_PATH=data/annotations/variant_summary.txt.gz
ALPHAMISSENSE_PATH=data/annotations/AlphaMissense_hg38.tsv.gz

# Model Settings
LLM_PROVIDER=anthropic
LLM_MODEL=claude-sonnet-4-20250514
EMBEDDING_MODEL=BAAI/bge-small-en-v1.5

# Performance
BATCH_SIZE=1000
TOP_K_RESULTS=20
```

### LLM Options

```bash
# Option 1: Anthropic Claude (Recommended)
LLM_PROVIDER=anthropic
LLM_MODEL=claude-sonnet-4-20250514

# Option 2: Local Ollama
LLM_PROVIDER=ollama
LLM_MODEL=llama3.1:70b

# Option 3: OpenAI
LLM_PROVIDER=openai
LLM_MODEL=gpt-4-turbo
```

---

## Directory Structure

```
rag-chat-pipeline/
├── docker-compose.yml          # Milvus + Attu services
├── requirements.txt            # Python dependencies
├── run.sh                      # Main CLI
├── .env                        # Configuration
│
├── config/
│   └── settings.py             # Application settings
│
├── src/
│   ├── vcf_parser.py           # VCF → Evidence objects
│   ├── annotator.py            # ClinVar + AlphaMissense annotation
│   ├── embedder.py             # Text → Vectors (BGE)
│   ├── milvus_client.py        # Vector DB operations
│   ├── llm_client.py           # LLM providers
│   ├── rag_engine.py           # RAG orchestration
│   ├── knowledge.py            # Clinker knowledge base (201 genes)
│   └── target_hypothesis.py    # Export to drug discovery
│
├── scripts/
│   ├── ingest_vcf.py           # Ingestion script
│   └── run_chat.py             # Start chat UI
│
├── app/
│   └── chat_ui.py              # Streamlit interface
│
└── data/
    ├── annotations/            # ClinVar, AlphaMissense
    │   ├── variant_summary.txt.gz
    │   └── AlphaMissense_hg38.tsv.gz
    ├── input/                  # VCF files (File Manager: browse/upload)
    │   └── HG002.genome.vcf.gz
    ├── output/                 # Results, exports (File Manager: browse/download)
    │   └── targets_for_phase5.json
    ├── targets/                # Saved target hypotheses
    └── cache/                  # Embedding cache
```

---

## Services

| Service | Port | Description |
|---------|------|-------------|
| **Streamlit** | 8501 | Chat interface + File Manager |
| **Milvus** | 19530 | Vector database |
| **Attu** | 8000 | Milvus web UI |

### UI Sidebar Navigation

The Streamlit interface includes a sidebar with multiple sections:

| Tab | Function |
|-----|----------|
| **Filters** | Search filters by gene, chromosome, impact level |
| **Targets** | View and manage target hypotheses |
| **Files** | File Manager - browse/upload VCF files |
| **VCF Preview** | Preview VCF file contents |
| **Metrics** | LLM performance metrics (TTFT, tokens/sec) |

---

## Troubleshooting

### Milvus Won't Start

```bash
docker-compose logs milvus
docker-compose down -v  # Reset volumes
docker-compose up -d
```

### Out of Memory During Ingestion

```bash
# Reduce batch size
python scripts/ingest_vcf.py --batch-size 100

# Or ingest annotated only
python scripts/ingest_vcf.py --annotated-only
```

### Anthropic API Error

```bash
# Verify API key
echo $ANTHROPIC_API_KEY

# Test connection
curl https://api.anthropic.com/v1/messages \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "content-type: application/json" \
  -d '{"model": "claude-sonnet-4-20250514", "max_tokens": 10, "messages": [{"role": "user", "content": "Hi"}]}'
```

### No Results for Query

- Ensure ingestion completed successfully
- Check if Milvus collection is loaded: `./run.sh status`
- Try more specific queries (gene names work better than vague descriptions)

---

## Related Pipelines

| Stage | Pipeline | Description |
|-------|----------|-------------|
| **1** | [Genomics Pipeline](https://github.com/ajones1923/hcls-ai-factory/tree/main/genomics-pipeline) | FASTQ → VCF with Parabricks |
| **2** | **RAG/Chat Pipeline** (This repo) | VCF → Target Hypothesis |
| **3** | [Drug Discovery Pipeline](https://github.com/ajones1923/hcls-ai-factory/tree/main/drug-discovery-pipeline) | Target → Molecule Candidates |

### Integration Flow

```
Genomics Pipeline          RAG/Chat Pipeline          Drug Discovery
     │                           │                          │
     │  HG002.genome.vcf.gz     │                          │
     └──────────────────────────▶│                          │
                                 │  Target Hypothesis       │
                                 │  (VCP, BRCA1, etc.)     │
                                 └──────────────────────────▶│
                                                            │
                                                   Molecule Candidates
```

---

## References

### Databases

- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [AlphaMissense](https://github.com/google-deepmind/alphamissense)
- [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/)

### Tools

- [Milvus Vector Database](https://milvus.io/)
- [BGE Embeddings](https://huggingface.co/BAAI/bge-small-en-v1.5)
- [Anthropic Claude](https://www.anthropic.com/)
- [Streamlit](https://streamlit.io/)

### Related Projects

- [NVIDIA DGX Spark](https://www.nvidia.com/en-us/data-center/dgx-spark/)
- [NVIDIA Clara](https://www.nvidia.com/en-us/clara/)
- [NVIDIA BioNeMo](https://www.nvidia.com/en-us/clara/bionemo/)

---

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](../LICENSE) file for details.

---

## Acknowledgments

- **NVIDIA** for DGX Spark and Clara ecosystem
- **Anthropic** for Claude API
- **Google DeepMind** for AlphaMissense predictions
- **NCBI** for ClinVar database
- **Ensembl** for VEP annotations
- **Milvus** for vector database technology

---

**Note:** This pipeline uses the GIAB HG002 reference genome for demonstration. For clinical use, ensure compliance with relevant regulations and validation requirements.
