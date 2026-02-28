# Precision Oncology Agent -- Architecture Design Document

**Author:** Adam Jones
**Date:** February 2026
**Version:** 1.0.0
**License:** Apache 2.0

---

## 1. Executive Summary

The Precision Oncology Agent extends the HCLS AI Factory platform to support closed-loop clinical decision support for molecular tumor boards. It transforms raw genomic data (VCF files from Stage 1 Parabricks) into structured treatment recommendations by combining:

1. **Variant Annotation** -- Classifying somatic and germline variants against ~40 actionable targets using AMP/ASCO/CAP evidence tiers
2. **Evidence Retrieval** -- Multi-collection RAG across 11 Milvus collections (~1,490+ owned vectors + 3.5M shared genomic vectors)
3. **Therapy Ranking** -- Evidence-level-sorted treatment recommendations with resistance awareness and contraindication flags
4. **Trial Matching** -- Hybrid deterministic + semantic clinical trial matching with composite scoring
5. **MTB Packet Generation** -- Structured Molecular Tumor Board packets exported as Markdown, JSON, PDF, or FHIR R4

### Key Results

| Metric | Value |
|---|---|
| Total vectors indexed | **~1,490+** across 10 owned collections + **3.5M** genomic vectors (read-only) |
| Multi-collection search latency | **< 200 ms** (11 collections, top-5 each) |
| Comparative dual retrieval | **~400 ms** (2 x 11 collections, entity-filtered) |
| Full RAG query (search + Claude) | **~24 sec** end-to-end |
| MTB packet generation | **< 30 sec** (VCF parse + annotation + evidence + ranking + trial match) |
| Knowledge graph entities | **~135+** (40 targets, 30 therapies, 20 resistance, 10 pathways, 15 biomarkers, 50+ aliases) |
| Query expansion | **12 maps**, ~120 keywords -> ~700 terms |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Precision Oncology Agent Role |
|---|---|
| **DataStore** | Raw files: VCF (Stage 1), PubMed XML, ClinicalTrials.gov JSON, CIViC JSON, seed data |
| **DataEngine** | Event-driven ingest pipelines (fetch -> parse -> embed -> store) for 8 data sources |
| **DataBase** | 11 Milvus collections (10 owned + 1 read-only) + knowledge graph (40 targets, 30 therapies, 20 resistance) |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + query expansion + case manager + trial matcher + therapy ranker |
| **AgentEngine** | OncoIntelligenceAgent (plan -> search -> evaluate -> synthesize) + Streamlit MTB Workbench + FastAPI REST |

### 2.2 System Diagram

```
                        +----------------------------------+
                        |  Streamlit MTB Workbench (8526)   |
                        |  5 tabs: Case | Evidence | Trial  |
                        |  | Therapy Ranker | Outcomes      |
                        +-----------------+----------------+
                                          |
                        +-----------------v----------------+
                        |     FastAPI REST Server (8527)    |
                        |  /query  /search  /api/cases      |
                        |  /api/trials/match  /api/reports   |
                        +-----------------+----------------+
                                          |
                        +-----------------v----------------+
                        |     OncoRAGEngine                 |
                        |  retrieve -> augment -> generate  |
                        |  + comparative detection          |
                        +-----------------+----------------+
                                          |
                  +--------------  "X vs Y"?  ---------------+
                  | YES                                  NO  |
                  v                                          v
        +------------------+                   +------------------+
        | Comparative Mode |                   | Standard Mode    |
        | Parse 2 entities |                   | Single retrieve  |
        | Dual retrieval   |                   |                  |
        +--------+---------+                   +--------+---------+
                 |                                      |
                 +-----------------+--------------------+
                                   |
                 +---------+-------+--------+---------+
                 |         |                |         |
        +--------v---+ +---v----------+ +--v------+ +v-----------+
        |  Query     | | Knowledge    | | Case    | | Claude     |
        |  Expansion | | Graph        | | Manager | | Sonnet 4.6 |
        |  12 maps   | | 40 targets   | | VCF->   | | Streaming  |
        |  ~700 terms| | 30 therapies | | MTB     | | RAG        |
        |            | | 20 resistance| | packets | |            |
        +--------+---+ +------+-------+ +---+----+ +------------+
                 |             |              |
        +--------v-------------v--------------v--------+
        |        Multi-Collection RAG Engine            |
        |   Parallel search across 11 collections      |
        |   Weighted: variant 0.18 | lit 0.16 |         |
        |   therapy 0.14 | guideline 0.12 |             |
        |   trial 0.10 | biomarker 0.08 | ...           |
        +--+----+-----+-----+------+------+---+-------+
           |    |     |     |      |      |   |
    +------v+ +v----++v----++v---+ +v----++v--v--+ +v---------+
    | onco_ | |onco_||onco_||onco|| onco_||onco_ | |genomic_  |
    | vari- | |lite-||ther-||guid|| trial||biom- | |evidence  |
    | ants  | |ratu-||apie-||elin||  s   ||arker | |(read-    |
    | ~300  | |re   ||s    ||es  || ~200 ||s ~80 | | only)    |
    |       | |~500 ||~120 ||~100||      ||      | | 3.5M     |
    +-------+ +-----++-----++----+ +-----++------+ +----------+
        ^        ^       ^      ^       ^      ^       ^
    +---+----+ +-+--+ +-+-+ +--+-+ +--+-+ +--+-+ +---+----+
    | CIViC  | |Pub | |FDA| |NCCN| |CT  | |TMB | |Stage 1 |
    | OncoKB | |Med | |   | |ASCO| |.gov| |MSI | |Parabr- |
    |        | |    | |   | |ESMO| |    | |HRD | |icks    |
    +--------+ +----+ +---+ +----+ +----+ +----+ +--------+
```

---

## 3. Data Collections

All 11 collections (10 owned + 1 read-only) use BGE-small-en-v1.5 embeddings (dim=384) with IVF_FLAT/COSINE indexing.

### 3.1 `onco_variants` -- ~300 records

| Attribute | Value |
|---|---|
| **Source** | CIViC database, OncoKB, curated seed data |
| **Fields** | id, gene, variant_name, variant_type, cancer_type, evidence_level, drugs, civic_id, vrs_id, text_summary, clinical_significance, allele_frequency |
| **Embedding** | FLOAT_VECTOR(384), BGE-small-en-v1.5 |
| **Index** | IVF_FLAT, COSINE, nlist=1024, nprobe=16 |
| **Evidence levels** | AMP/ASCO/CAP Tier A through D |
| **Variant types** | SNV, indel, CNV amplification, CNV deletion, fusion, rearrangement, structural variant |

### 3.2 `onco_literature` -- ~500 records

| Attribute | Value |
|---|---|
| **Source** | PubMed via NCBI E-utilities (esearch + efetch) |
| **Fields** | id (PMID), title, text_chunk, source_type, year, cancer_type, gene, variant, keywords, journal |
| **Cancer type tagging** | Automated classification across 20+ cancer types |
| **Gene extraction** | Automated from title and abstract |

### 3.3 `onco_therapies` -- ~120 records

| Attribute | Value |
|---|---|
| **Source** | FDA-approved labels, curated seed data |
| **Fields** | id, drug_name, category, targets, approved_indications, resistance_mechanisms, evidence_level, text_summary, mechanism_of_action |
| **Categories** | Targeted, immunotherapy, chemotherapy, hormonal, combination, radiotherapy, cell therapy |

### 3.4 `onco_guidelines` -- ~100 records

| Attribute | Value |
|---|---|
| **Source** | NCCN, ASCO, ESMO, WHO, CAP/AMP guideline excerpts |
| **Fields** | id, org, cancer_type, version, year, key_recommendations, text_summary, evidence_level |
| **Organizations** | NCCN, ESMO, ASCO, WHO, CAP/AMP |

### 3.5 `onco_trials` -- ~200 records

| Attribute | Value |
|---|---|
| **Source** | ClinicalTrials.gov REST API v2 |
| **Fields** | id (NCT ID), title, text_summary, phase, status, sponsor, cancer_types, biomarker_criteria, enrollment, start_year, outcome_summary |
| **Phase distribution** | Early Phase 1 through Phase 4 |
| **Status** | Recruiting, completed, active, terminated, withdrawn |

### 3.6 `onco_biomarkers` -- ~80 records

| Attribute | Value |
|---|---|
| **Source** | Curated biomarker panels, FDA companion diagnostics |
| **Fields** | id, name, biomarker_type, cancer_types, predictive_value, testing_method, clinical_cutoff, text_summary, evidence_level |
| **Types** | Predictive, prognostic, diagnostic, monitoring, resistance, pharmacodynamic |

### 3.7 `onco_resistance` -- ~80 records

| Attribute | Value |
|---|---|
| **Source** | Published resistance literature, curated mechanisms |
| **Fields** | id, primary_therapy, gene, mechanism, bypass_pathway, alternative_therapies, text_summary |
| **Coverage** | EGFR TKI resistance (T790M, C797S, MET bypass), BRAF resistance (MAPK reactivation), immunotherapy resistance (antigen loss, beta-2-microglobulin), PARP inhibitor resistance (BRCA reversion), endocrine resistance (ESR1 mutations) |

### 3.8 `onco_pathways` -- ~50 records

| Attribute | Value |
|---|---|
| **Source** | Curated pathway data |
| **Fields** | id, name, key_genes, therapeutic_targets, cross_talk, text_summary |
| **Pathways** | MAPK, PI3K/AKT/mTOR, DDR, cell cycle, apoptosis, Wnt, Notch, Hedgehog, JAK/STAT, angiogenesis |

### 3.9 `onco_outcomes` -- ~50 records

| Attribute | Value |
|---|---|
| **Source** | Synthetic treatment outcomes for closed-loop learning |
| **Fields** | id, case_id, therapy, cancer_type, response (RECIST), duration_months, toxicities, biomarkers_at_baseline, text_summary |
| **Response categories** | Complete response, partial response, stable disease, progressive disease, not evaluable |

### 3.10 `onco_cases` -- ~10 records

| Attribute | Value |
|---|---|
| **Source** | Synthetic de-identified patient snapshots |
| **Fields** | id, patient_id, cancer_type, stage, variants, biomarkers, prior_therapies, text_summary |

### 3.11 `genomic_evidence` (read-only) -- 3,561,170 records

| Attribute | Value |
|---|---|
| **Source** | Stage 1 Parabricks pipeline (VCF -> annotation -> Milvus) |
| **Fields** | id, chrom, pos, ref, alt, qual, gene, consequence, impact, genotype, text_summary, clinical_significance, rsid, disease_associations, am_pathogenicity, am_class |
| **Shared** | Read-only -- populated by Stage 1, consumed by all agents |

### 3.12 Index Configuration (all collections)

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 1024 |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Knowledge Graph Detail

### 4.1 Actionable Targets (~40 entries)

Each entry includes: gene name, full name, cancer types, key variants, targeted therapies, combination therapies, resistance mutations, pathway, evidence level, and description.

| Target | Cancer Types | Key Therapies |
|---|---|---|
| **BRAF** | Melanoma, NSCLC, CRC, thyroid | vemurafenib, dabrafenib+trametinib, encorafenib+binimetinib |
| **EGFR** | NSCLC | osimertinib, erlotinib, gefitinib, amivantamab |
| **ALK** | NSCLC | alectinib, lorlatinib, brigatinib, crizotinib |
| **KRAS** | NSCLC, CRC, pancreatic | sotorasib, adagrasib |
| **HER2/ERBB2** | Breast, gastric, NSCLC | trastuzumab, T-DXd, tucatinib |
| **ROS1** | NSCLC | crizotinib, entrectinib, lorlatinib |
| **NTRK** | Tissue-agnostic | larotrectinib, entrectinib |
| **RET** | NSCLC, thyroid, MTC | selpercatinib, pralsetinib |
| **MET** | NSCLC | capmatinib, tepotinib |
| **FGFR** | Urothelial, cholangiocarcinoma | erdafitinib, futibatinib, pemigatinib |
| **PIK3CA** | Breast | alpelisib |
| **BRCA1/2** | Breast, ovarian, prostate, pancreatic | olaparib, rucaparib, niraparib |
| **IDH1/2** | AML, cholangiocarcinoma, glioma | ivosidenib, enasidenib |
| + ~26 more | Various | Various |

### 4.2 Therapy Mapping (~30 entries)

| Drug | Brand | Category | Targets | Key Trials |
|---|---|---|---|---|
| osimertinib | Tagrisso | EGFR TKI (3rd gen) | EGFR | FLAURA, FLAURA2 |
| sotorasib | Lumakras | KRAS G12C inhibitor | KRAS G12C | CodeBreaK 200 |
| pembrolizumab | Keytruda | Anti-PD-1 | PD-L1, MSI-H, TMB-H | KEYNOTE-024/158/177 |
| dabrafenib+trametinib | Tafinlar+Mekinist | BRAF+MEK combo | BRAF V600 | COMBI-d, COMBI-v |
| lorlatinib | Lorbrena | ALK TKI (3rd gen) | ALK | CROWN |
| olaparib | Lynparza | PARP inhibitor | BRCA, HRD | SOLO-1, PROfound |
| trastuzumab deruxtecan | Enhertu | ADC | HER2 | DESTINY-Breast |
| + ~23 more | Various | Various | Various | Various |

### 4.3 Resistance Map (~20 entries)

| Target Drug | Resistance Mechanism | Bypass Pathway | Alternatives |
|---|---|---|---|
| Erlotinib/Gefitinib | EGFR T790M gatekeeper | -- | Osimertinib |
| Osimertinib | EGFR C797S, MET amplification | MET bypass | Amivantamab, MET inhibitor combo |
| Vemurafenib/Dabrafenib | NRAS activation, BRAF amplification | MAPK reactivation | MEK inhibitor combo, immunotherapy |
| Crizotinib (ALK) | ALK G1202R solvent-front | -- | Lorlatinib |
| Pembrolizumab | Beta-2-microglobulin loss, JAK1/2 mutation | Antigen presentation loss | Combination IO, CTLA-4 |
| Olaparib | BRCA reversion mutation, 53BP1 loss | HR restoration | Platinum re-challenge |
| Trastuzumab | HER2 truncation (p95-HER2), PI3K activation | HER3 upregulation | T-DXd, tucatinib |
| + ~13 more | Various | Various | Various |

### 4.4 Pathway Map (~10 entries)

| Pathway | Key Genes | Druggable Nodes | Cross-Talk |
|---|---|---|---|
| MAPK | KRAS, NRAS, BRAF, MEK1/2, ERK | BRAF, MEK, ERK | PI3K, cell cycle |
| PI3K/AKT/mTOR | PIK3CA, PTEN, AKT1, mTOR | PI3K, AKT, mTOR | MAPK, apoptosis |
| DDR | BRCA1, BRCA2, ATM, ATR, PALB2 | PARP, ATR, CHK1 | Cell cycle, apoptosis |
| Cell Cycle | CDK4, CDK6, RB1, CDKN2A, CCND1 | CDK4/6 | MAPK, PI3K |
| Apoptosis | BCL2, BAX, BAK, MCL1 | BCL-2, MCL-1 | DDR, cell cycle |
| Angiogenesis | VEGFA, VEGFR, PDGFR, FGF | VEGFR, PDGFR | MAPK, Notch |
| + 4 more | Various | Various | Various |

### 4.5 Biomarker Panels (~15 entries)

| Biomarker | Type | Cutoff | Therapy Implication |
|---|---|---|---|
| TMB-H | Predictive | >= 10 mut/Mb | Pembrolizumab (KEYNOTE-158) |
| MSI-H / dMMR | Predictive | IHC / PCR / NGS | Pembrolizumab (tissue-agnostic) |
| PD-L1 TPS | Predictive | >= 50% (NSCLC 1L) | Pembrolizumab monotherapy |
| PD-L1 CPS | Predictive | >= 10 (gastric) | Pembrolizumab + chemo |
| HRD | Predictive | GIS score | PARP inhibitors |
| ALK rearrangement | Predictive | FISH / IHC / NGS | ALK TKIs |
| NTRK fusion | Predictive | IHC / FISH / NGS | Larotrectinib, entrectinib |
| ROS1 fusion | Predictive | FISH / NGS | Crizotinib, entrectinib |
| EGFR mutation | Predictive | NGS / PCR | EGFR TKIs |
| BRCA mutation | Predictive | NGS (germline + somatic) | PARP inhibitors |
| ctDNA / cfDNA | Monitoring | VAF dynamics | MRD tracking, resistance detection |
| FGFR alteration | Predictive | NGS | FGFR inhibitors |
| RET fusion | Predictive | NGS / FISH | Selpercatinib, pralsetinib |
| BRAF V600E | Predictive | NGS / IHC | BRAF+MEK combo |
| KRAS G12C | Predictive | NGS | Sotorasib, adagrasib |

---

## 5. Query Expansion

Twelve expansion map categories for improved RAG retrieval recall:

| Category | Keywords | Expanded Terms | Examples |
|---|---|---|---|
| Cancer Type | 10 | ~50 | NSCLC -> [lung adenocarcinoma, lung squamous, EGFR-mutant lung, ALK-positive lung, ...] |
| Gene | 10 | ~50 | EGFR -> [epidermal growth factor receptor, EGFR L858R, EGFR exon 19 deletion, T790M, C797S, ...] |
| Therapy | 10 | ~50 | osimertinib -> [Tagrisso, EGFR TKI, third-generation EGFR inhibitor, FLAURA, FLAURA2, ...] |
| Biomarker | 10 | ~55 | MSI-H -> [microsatellite instability high, mismatch repair deficient, dMMR, MLH1 loss, ...] |
| Pathway | 10 | ~50 | MAPK -> [RAS-RAF-MEK-ERK, MAPK cascade, ERK signaling, RAS signaling, ...] |
| Resistance | 10 | ~55 | EGFR resistance -> [T790M gatekeeper, C797S, MET amplification bypass, small cell transformation, ...] |
| Clinical | 10 | ~50 | response -> [overall response rate, ORR, complete response, partial response, RECIST criteria, ...] |
| Trial | 10 | ~50 | basket trial -> [histology-agnostic, biomarker-selected, tissue-agnostic, master protocol, ...] |
| Immunotherapy | 10 | ~50 | immunotherapy -> [checkpoint inhibitor, anti-PD-1, anti-PD-L1, anti-CTLA-4, ICB, ...] |
| Surgery/Radiation | 10 | ~50 | neoadjuvant -> [preoperative, downstaging, pathologic complete response, pCR, ...] |
| Toxicity | 10 | ~50 | irAE -> [immune-related adverse event, checkpoint toxicity, autoimmune toxicity, steroid taper, ...] |
| Genomics | 10 | ~50 | sequencing -> [next-generation sequencing, NGS, whole exome, targeted panel, liquid biopsy, ctDNA, ...] |
| **Total** | **~120** | **~700** | |

The `expand_query()` function detects keywords in the user's query (case-insensitive substring match) and returns up to 10 semantically related expansion terms, which are used to run additional filtered searches across collections.

---

## 6. RAG Engine

### 6.1 Collection Weights

| Collection | Weight | Rationale |
|---|---|---|
| onco_variants | 0.18 | Actionable variant evidence is the foundation of precision oncology |
| onco_literature | 0.16 | Published evidence provides clinical context and validation |
| onco_therapies | 0.14 | Therapy data directly informs treatment recommendations |
| onco_guidelines | 0.12 | NCCN/ASCO/ESMO guidelines represent expert consensus |
| onco_trials | 0.10 | Active trials offer options when standard therapy fails |
| onco_biomarkers | 0.08 | Biomarker status drives therapy selection (TMB, MSI, PD-L1) |
| onco_resistance | 0.07 | Resistance awareness prevents ineffective re-treatment |
| onco_pathways | 0.06 | Pathway context explains mechanism and cross-talk |
| onco_outcomes | 0.04 | Historical outcomes calibrate expectations |
| onco_cases | 0.02 | Similar cases provide clinical analogy |
| genomic_evidence | 0.03 | Patient-specific variant context from Stage 1 |

### 6.2 Retrieval Flow

```
User Query: "What therapies target BRAF V600E in melanoma?"
    |
    +-- 1. Embed query (BGE asymmetric prefix)                         [< 5 ms]
    |      "Represent this sentence for searching relevant passages: ..."
    |
    +-- 2. Parallel search across 11 collections (top-5 each)         [< 200 ms]
    |   +-- onco_variants:     BRAF V600E records                      (score: 0.85-0.92)
    |   +-- onco_literature:   BRAF melanoma papers                    (score: 0.78-0.88)
    |   +-- onco_therapies:    vemurafenib, dabrafenib, encorafenib    (score: 0.82-0.90)
    |   +-- onco_guidelines:   NCCN melanoma guidelines                (score: 0.80-0.87)
    |   +-- onco_trials:       BRAF V600E melanoma trials              (score: 0.75-0.85)
    |   +-- onco_resistance:   BRAF resistance mechanisms              (score: 0.72-0.82)
    |
    +-- 3. Query expansion: "BRAF" -> [B-Raf proto-oncogene, V600E,    [< 1 ms]
    |      V600K, class I mutation, BRAF fusion, ...]
    |
    +-- 4. Expanded filtered search (top-3 per expanded term)          [< 200 ms]
    |
    +-- 5. Merge + deduplicate + weighted rank (cap at 30 results)     [< 1 ms]
    |
    +-- 6. Knowledge graph augmentation:                               [< 1 ms]
    |      BRAF -> targeted_therapies: [vemurafenib, dabrafenib, encorafenib]
    |      BRAF -> combination: [dabrafenib+trametinib, encorafenib+binimetinib]
    |      BRAF -> resistance: [NRAS activation, BRAF amplification, MAP2K1 C121S]
    |      BRAF -> pathway: MAPK
    |
    +-- 7. Build prompt: evidence grouped by collection +              [< 1 ms]
    |      knowledge context + question + citation instructions
    |
    +-- 8. Stream Claude Sonnet 4.6 response                           [~22-24 sec]
           Grounded answer with clickable PubMed and ClinicalTrials.gov citations
```

### 6.3 Comparative Analysis

Comparative queries auto-detected via regex: `compare`, `vs`, `versus`, `difference between`, `head.to.head`.

```
User: "Compare EGFR TKI generations for NSCLC"
    |
    +-- 1. _is_comparative() detects "compare"                         [< 1 ms]
    +-- 2. _parse_comparison_entities() -> "EGFR TKI generations", "NSCLC"
    +-- 3. Dual retrieve() for each entity                             [~400 ms]
    +-- 4. Identify shared / head-to-head evidence                     [< 1 ms]
    +-- 5. _build_comparative_prompt() with structured sections        [< 1 ms]
           - Evidence for Entity A, Evidence for Entity B
           - Shared / Head-to-Head Evidence
           - Instructions: MoA, efficacy, safety, biomarkers,
             resistance, guidelines, trials, recommendation
    +-- 6. Stream Claude Sonnet 4.6 (max_tokens=4096)                  [~28-30 sec]
           Structured comparison with tables, pros/cons, clinical context
```

---

## 7. Case Management Workflow

### OncologyCaseManager

The case manager handles the complete patient case lifecycle:

```
Step 1: VCF Input
    +-- Accept raw VCF text or pre-parsed variant list
    +-- Parse: extract CHROM, POS, REF, ALT, GENE, consequence from INFO field
    +-- Filter: PASS variants only

Step 2: Variant Annotation
    +-- For each variant, classify actionability against ACTIONABLE_TARGETS
    +-- AMP/ASCO/CAP evidence tiers:
         A = FDA-approved companion diagnostic
         B = Clinical evidence from well-powered studies
         C = Case reports and small series
         D = Preclinical / in-vitro data
    +-- Map to drugs, resistance mechanisms, pathway context

Step 3: Case Snapshot
    +-- Create CaseSnapshot with: patient_id, cancer_type, stage,
        variants (with actionability), biomarkers, prior_therapies
    +-- Embed case summary and store in onco_cases collection

Step 4: MTB Packet Generation
    +-- Build variant table (all variants with evidence levels)
    +-- Build evidence table (RAG retrieval per actionable variant)
    +-- Build therapy ranking (evidence-level sort + resistance + contraindication)
    +-- Build trial matches (hybrid deterministic + semantic search)
    +-- Build open questions (VUS, missing biomarkers, evidence gaps)
    +-- Export as Markdown, JSON, PDF, or FHIR R4
```

---

## 8. Trial Matching Algorithm

The `TrialMatcher` uses a 4-step hybrid approach:

### Step 1: Deterministic Filter

Filter `onco_trials` by `cancer_type == patient.cancer_type AND status == "Recruiting"`. Returns up to 30 deterministic hits.

### Step 2: Semantic Search

Build eligibility query from patient profile: `"{cancer_type} clinical trial stage {stage} {biomarkers}"`. Embed and search `onco_trials` for top 30 semantic hits.

### Step 3: Composite Scoring

For each trial, compute:

```
composite_score = 0.40 * biomarker_match_score
                + 0.25 * semantic_score
                + 0.20 * phase_weight
                + 0.15 * status_weight
```

| Component | Range | Description |
|---|---|---|
| biomarker_match_score | 0.0 - 1.0 | Fraction of patient biomarkers mentioned in trial criteria |
| semantic_score | 0.0 - 1.0 | Cosine similarity from vector search |
| phase_weight | Phase 3=1.0, Phase 2/3=0.9, Phase 2=0.8, Phase 1/2=0.7, Phase 1=0.6 | Higher phase = more relevant |
| status_weight | Recruiting=1.0, Enrolling by invitation=0.8, Active not recruiting=0.6, Not yet=0.4 | Recruiting preferred |

### Step 4: Explanation Generation

For each matched trial, generate structured explanation:
- Matched criteria (cancer type, biomarkers found in trial criteria)
- Unmatched criteria (biomarkers not confirmed)
- Trial phase and status summary

---

## 9. Therapy Ranking System

The `TherapyRanker` implements a 6-step evidence-based ranking:

### Step 1: Variant-Driven Therapy Identification

For each variant, check `ACTIONABLE_TARGETS` for matching gene/variant. Extract targeted therapies with evidence levels.

### Step 2: Biomarker-Driven Therapy Identification

Check standard biomarker-therapy mappings:

| Biomarker | Threshold | Therapy | Evidence |
|---|---|---|---|
| MSI-H / dMMR | Positive | Pembrolizumab | Level A |
| TMB-H | >= 10 mut/Mb | Pembrolizumab | Level A |
| HRD / BRCA | Positive | PARP inhibitors | Level A/B |
| PD-L1 TPS | >= 50% | Pembrolizumab 1L | Level A |
| NTRK fusion | Positive | Larotrectinib, entrectinib | Level A |

### Step 3: Evidence-Level Sorting

Rank all candidate therapies by AMP/ASCO/CAP evidence level: A > B > C > D > E.

### Step 4: Resistance Check

Cross-reference candidate drugs against `RESISTANCE_MAP` and patient's prior therapies. Flag drugs where prior treatment may have induced resistance.

### Step 5: Contraindication Check

Flag drugs in the same drug class as a previously failed therapy (via `THERAPY_MAP` drug class lookup).

### Step 6: Supporting Evidence Retrieval

For each therapy, search `onco_therapies` and `onco_literature` for supporting evidence. Attach top citations.

Final ranking: clean therapies first (sorted by evidence level), resistance/contraindicated therapies demoted to bottom of their tier.

---

## 10. Cross-Modal Integration

### Variant Severity -> Imaging Collections

When actionable variants with evidence level A or B are identified, the `OncoCrossModalTrigger` module:

1. Queries `genomic_evidence` collection (3.5M vectors) for variant context
2. Discovers and queries `imaging_*` collections (if Imaging Agent is deployed)
3. Correlates genomic findings with imaging characteristics (e.g., EGFR mutation + lung nodule imaging)

### Variant Actionability -> Drug Discovery

When the Precision Oncology Agent identifies an actionable target:

1. Target protein structure is identified (e.g., BRAF V600E -> PDB structure)
2. Cross-modal trigger passes target to Stage 3 Drug Discovery pipeline
3. BioNeMo MolMIM generates novel inhibitor candidates
4. DiffDock simulates binding, RDKit scores drug-likeness

### Graceful Degradation

If imaging or drug discovery collections are unavailable, cross-modal queries silently return empty results. The core oncology workflow continues unimpaired.

---

## 11. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 11.1 Search Performance

| Operation | Latency | Notes |
|---|---|---|
| Single collection search (top-5) | 3-8 ms | Milvus IVF_FLAT with cached index |
| 11-collection parallel search (top-5 each) | < 200 ms | ThreadPoolExecutor, up to 8 threads |
| Query expansion + filtered re-search | < 200 ms | Up to 10 expanded terms |
| Knowledge graph augmentation | < 1 ms | In-memory dictionary lookup |
| Full retrieve() pipeline | < 500 ms | Embed + search + expand + merge + knowledge |
| Comparative dual retrieval | ~400 ms | Two retrieve() calls |

### 11.2 RAG Query Performance

| Operation | Latency | Notes |
|---|---|---|
| Full query() (retrieve + Claude generate) | ~24 sec | Dominated by LLM generation |
| Comparative query (dual retrieve + Claude) | ~30 sec | Structured comparison prompt |
| Streaming query_stream() (time to first token) | ~3 sec | Evidence returned immediately |
| Response length (standard) | 800-2000 chars | Grounded answer with citations |
| Response length (comparative) | 1500-3000 chars | Structured tables + clinical context |

### 11.3 Case Management Performance

| Operation | Latency | Notes |
|---|---|---|
| VCF parsing (1000 variants) | < 500 ms | Regex extraction, PASS filter |
| Variant actionability classification | < 100 ms | In-memory knowledge graph lookup |
| MTB packet generation (full) | < 30 sec | VCF + annotation + evidence + ranking + trials |
| Trial matching (10 results) | < 10 sec | Deterministic + semantic + scoring |
| Therapy ranking (all candidates) | < 5 sec | Variant + biomarker + resistance + evidence |

---

## 12. Infrastructure

### 12.1 Technology Stack

| Component | Technology | Version/Detail |
|---|---|---|
| Vector database | Milvus | 2.4, localhost:19530 |
| Embedding model | BGE-small-en-v1.5 | 384-dim, BAAI, ~33M params |
| LLM | Claude Sonnet 4.6 | Anthropic API, claude-sonnet-4-20250514 |
| UI framework | Streamlit | Port 8526, 5-tab MTB Workbench |
| API framework | FastAPI | Port 8527, OpenAPI docs at /docs |
| Data models | Pydantic | BaseModel + Field validation |
| Configuration | Pydantic BaseSettings | ONCO_ prefix environment variables |
| PDF export | ReportLab | NVIDIA-themed Platypus layout |
| FHIR export | Custom | R4 Bundle with SNOMED CT + LOINC coding |
| Hardware target | NVIDIA DGX Spark | GB10 GPU, 128GB unified, $3,999 |

### 12.2 Docker Services

| Service | Image | Port | Role |
|---|---|---|---|
| `milvus-etcd` | quay.io/coreos/etcd:v3.5.5 | 2379 | Milvus metadata store |
| `milvus-minio` | minio/minio:v2023.03 | 9000, 9001 | Milvus object storage |
| `milvus-standalone` | milvusdb/milvus:v2.4 | 19530, 9091 | Vector database |
| `onco-streamlit` | Built from Dockerfile | 8526 | Streamlit MTB Workbench |
| `onco-api` | Built from Dockerfile | 8527 | FastAPI REST server |
| `onco-setup` | Built from Dockerfile | -- | One-shot collection setup + seed |

### 12.3 Service Ports

| Port | Service |
|---|---|
| 8526 | Precision Oncology MTB Workbench (Streamlit) |
| 8527 | Precision Oncology REST API (FastAPI) |
| 19530 | Milvus vector database (shared with main pipeline) |

---

## 13. File Structure

```
precision_oncology_agent/agent/
+-- src/
|   +-- __init__.py
|   +-- models.py                          # Pydantic data models + MTBPacket (497 lines)
|   +-- collections.py                     # Milvus collection schemas + manager (606 lines)
|   +-- knowledge.py                       # Knowledge graph: targets, therapies, resistance (1,194 lines)
|   +-- query_expansion.py                 # 12 expansion maps, ~120 -> ~700 terms (676 lines)
|   +-- rag_engine.py                      # Multi-collection RAG + comparative (780 lines)
|   +-- agent.py                           # Plan-search-synthesize pipeline (489 lines)
|   +-- case_manager.py                    # VCF parsing + MTB packet generation (509 lines)
|   +-- trial_matcher.py                   # Hybrid trial matching (393 lines)
|   +-- therapy_ranker.py                  # Evidence-based therapy ranking (552 lines)
|   +-- cross_modal.py                     # Cross-modal triggers (395 lines)
|   +-- export.py                          # MD, JSON, PDF, FHIR R4 export (876 lines)
|   +-- metrics.py                         # Prometheus metrics (362 lines)
|   +-- scheduler.py                       # Data ingestion scheduler (263 lines)
|   +-- ingest/
|   |   +-- __init__.py
|   |   +-- base.py                        # Base ingest pipeline (249 lines)
|   |   +-- civic_parser.py               # CIViC variant ingest (340 lines)
|   |   +-- oncokb_parser.py              # OncoKB parser (104 lines)
|   |   +-- literature_parser.py          # PubMed E-utilities ingest (248 lines)
|   |   +-- clinical_trials_parser.py     # ClinicalTrials.gov API v2 (279 lines)
|   |   +-- guideline_parser.py           # NCCN/ASCO/ESMO parser (168 lines)
|   |   +-- pathway_parser.py             # Pathway parser (121 lines)
|   |   +-- resistance_parser.py          # Resistance parser (125 lines)
|   |   +-- outcome_parser.py             # Outcome parser (158 lines)
|   +-- utils/
|       +-- __init__.py
|       +-- vcf_parser.py                 # VCF parsing utilities (361 lines)
|       +-- pubmed_client.py              # NCBI E-utilities HTTP client (296 lines)
+-- app/
|   +-- oncology_ui.py                    # Streamlit MTB Workbench (703 lines)
+-- api/
|   +-- __init__.py
|   +-- main.py                           # FastAPI application (347 lines)
|   +-- routes/
|       +-- __init__.py
|       +-- meta_agent.py                 # /api/ask, /api/deep-research (169 lines)
|       +-- cases.py                      # /api/cases endpoints (234 lines)
|       +-- trials.py                     # /api/trials/match (153 lines)
|       +-- reports.py                    # /api/reports/{format} (236 lines)
|       +-- events.py                     # /api/events, cross-modal (89 lines)
+-- config/
|   +-- settings.py                       # Pydantic BaseSettings (109 lines)
+-- tests/
|   +-- __init__.py
|   +-- conftest.py                       # Test fixtures (214 lines)
+-- requirements.txt
+-- LICENSE                               # Apache 2.0
```

**39 Python files | ~12,301 lines of code | Apache 2.0**

---

## 14. Implementation Status

| File | Lines | Status |
|---|---|---|
| src/models.py | 497 | Complete -- 12 enums, 10 domain models, MTBPacket, search models |
| src/collections.py | 606 | Complete -- 11 collection schemas, OncoCollectionManager |
| src/knowledge.py | 1,194 | Complete -- ~40 targets, ~30 therapies, ~20 resistance, ~10 pathways, ~15 biomarkers |
| src/query_expansion.py | 676 | Complete -- 12 expansion maps, ~120 keywords -> ~700 terms |
| src/rag_engine.py | 780 | Complete -- OncoRAGEngine, comparative mode, citation formatting |
| src/agent.py | 489 | Complete -- OncoIntelligenceAgent, plan-search-evaluate-synthesize |
| src/case_manager.py | 509 | Complete -- VCF parsing, case lifecycle, MTB packet generation |
| src/trial_matcher.py | 393 | Complete -- 4-step hybrid matching with composite scoring |
| src/therapy_ranker.py | 552 | Complete -- 6-step evidence-based ranking with resistance |
| src/cross_modal.py | 395 | Complete -- Genomic + imaging cross-modal triggers |
| src/export.py | 876 | Complete -- Markdown, JSON, PDF (ReportLab), FHIR R4 |
| src/metrics.py | 362 | Complete -- Prometheus metrics |
| src/scheduler.py | 263 | Complete -- APScheduler data refresh |
| src/ingest/*.py | 1,792 | Complete -- 8 ingest pipelines (base + 7 domain parsers) |
| src/utils/*.py | 657 | Complete -- VCF parser, PubMed client |
| app/oncology_ui.py | 703 | Complete -- 5-tab Streamlit MTB Workbench |
| api/main.py | 347 | Complete -- FastAPI app with 7 core endpoints |
| api/routes/*.py | 881 | Complete -- 5 route modules |
| config/settings.py | 109 | Complete -- OncoSettings with ONCO_ prefix |
| tests/conftest.py | 214 | Complete -- Test fixtures |

---

## 15. Credits

- **Adam Jones** -- HCLS AI Factory platform, 14+ years genomic research
- **Apache 2.0 License**
