# CAR-T Intelligence Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** February 2026
**Version:** 1.2.0
**License:** Apache 2.0

---

## 1. Executive Summary

The CAR-T Intelligence Agent extends the HCLS AI Factory platform to support cross-functional intelligence across the CAR-T cell therapy development lifecycle. The agent breaks down data silos between the stages of CAR-T development:

1. **Target Identification** — Antigen biology, expression profiling, disease association
2. **CAR Design** — scFv selection, costimulatory domains, signaling architecture
3. **Vector Engineering** — Transduction, viral vector production, manufacturing processes
4. **In Vitro / In Vivo Testing** — Cytotoxicity, cytokine assays, animal models, persistence
5. **Clinical Development** — Trial design, response rates, toxicity management

The platform enables cross-functional queries like *"Why do CD19 CAR-T therapies fail in relapsed B-ALL?"* that simultaneously search published literature, clinical trials, CAR construct data, assay results, and manufacturing records — returning grounded answers with clickable [PubMed](https://pubmed.ncbi.nlm.nih.gov/) and [ClinicalTrials.gov](https://clinicaltrials.gov/) citations.

**Comparative Analysis Mode** auto-detects "X vs Y" queries (e.g., *"Compare 4-1BB vs CD28 costimulatory domains"*), runs dual retrievals with per-entity filtering, and produces structured side-by-side analysis with markdown tables.

### Key Results

| Metric | Value |
|---|---|
| Total vectors indexed | **6,266** across 11 Milvus collections (10 owned + 1 read-only) |
| Multi-collection search latency | **12-16 ms** (11 collections, top-5 each, cached) |
| Comparative dual retrieval | **~365 ms** (2 × 11 collections, entity-filtered) |
| Full RAG query (search + Claude) | **~24 sec** end-to-end |
| Comparative RAG query (dual search + Claude) | **~30 sec** end-to-end |
| Cosine similarity scores | **0.74 - 0.90** on demo queries |
| Manufacturing success rate (seed script) | **100%** (all collections populated, 0 ingest errors) |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | CAR-T Agent Role |
|---|---|
| **DataStore** | Raw files: PubMed XML, ClinicalTrials.gov JSON, seed data JSON |
| **DataEngine** | Event-driven ingest pipelines (fetch → parse → embed → store) |
| **DataBase** | 11 Milvus collections (10 owned + 1 read-only) + knowledge graph (25 targets, 8 toxicities, 10 mfg) |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + query expansion |
| **AgentEngine** | CARTRAGEngine (retrieve → augment → generate) + Streamlit UI |

### 2.3 System Diagram

```
                        ┌─────────────────────────────┐
                        │   Streamlit Chat UI (8521)   │
                        │   Cross-functional queries   │
                        │   + Comparative Analysis UI  │
                        └──────────────┬──────────────┘
                                       │
                        ┌──────────────▼──────────────┐
                        │     CARTRAGEngine            │
                        │  retrieve → augment → gen    │
                        │  + comparative detection     │
                        └──────────────┬──────────────┘
                                       │
                  ┌────────────── "X vs Y"? ──────────────┐
                  │ YES                                NO  │
                  ▼                                        ▼
        ┌──────────────────┐                   ┌──────────────────┐
        │ Comparative Mode │                   │ Standard Mode    │
        │ Parse 2 entities │                   │ Single retrieve  │
        │ Entity resolution│                   │                  │
        │ Dual retrieval   │                   │                  │
        └────────┬─────────┘                   └────────┬─────────┘
                 │                                      │
                 └──────────────────┬───────────────────┘
                                   │
                 ┌─────────────────┼─────────────────────┐
                 │                 │                      │
        ┌────────▼────────┐  ┌────▼───────────┐  ┌──────▼──────────┐
        │  Query Expansion │  │ Knowledge Graph │  │ Claude Sonnet   │
        │  169 keywords    │  │ 25 targets      │  │ 4.6 (Anthropic) │
        │  → 1,496 terms   │  │ 8 toxicities    │  │ Streaming RAG   │
        │  12 categories   │  │ 10 manufacturing │  │ + Comparative   │
        │                  │  │ 39+ entity alias │  │   prompt builder│
        └────────┬────────┘  └────────┬────────┘  └─────────────────┘
                 │                     │
        ┌────────▼─────────────────────▼────────┐
        │        Multi-Collection RAG Engine     │
        │   Parallel search across 11 collections│
        │   Weighted: lit 0.30 | trial 0.25 |    │
        │   construct 0.20 | assay 0.15 |        │
        │   manufacturing 0.10                   │
        └───┬────┬─────┬─────┬──────┬───────────┘
            │    │     │     │      │
    ┌───────▼┐ ┌▼────┐┌▼────┐┌▼───┐┌▼────────┐
    │ cart_  │ │cart_ ││cart_││cart_││ cart_    │
    │ litera-│ │trial-││cons-││assa││ manufac- │
    │ ture   │ │s     ││truc-││ys  ││ turing   │
    │ 5,047  │ │ 973  ││ts 6 ││ 45 ││   30     │
    └────────┘ └─────┘└─────┘└────┘└──────────┘
         ▲        ▲       ▲      ▲       ▲
    ┌────┴────┐ ┌─┴──┐ ┌──┴──┐┌─┴──┐ ┌──┴───┐
    │ PubMed  │ │ CT │ │ FDA │ │Pub │ │ Pub  │
    │E-utils  │ │.gov│ │seed │ │lit │ │ lit  │
    │ API     │ │API │ │data │ │seed│ │ seed │
    └─────────┘ └────┘ └─────┘└────┘ └──────┘
```

---

## 3. Data Collections — Actual State

All 11 collections (10 owned + 1 read-only) are populated and searchable.

### 3.1 `cart_literature` — 5,047 records

| Attribute | Value |
|---|---|
| **Source** | PubMed via NCBI E-utilities (esearch + efetch) |
| **Ingest time** | ~15 min |
| **Fields** | PMID, title, text_chunk, source_type, year, cart_stage, target_antigen, disease, keywords, journal |
| **Embedding** | FLOAT_VECTOR(384), BGE-small-en-v1.5 |
| **Index** | IVF_FLAT, COSINE, nlist=1024, nprobe=16 |
| **Stage classification** | Automated: target_id, car_design, vector_eng, testing, clinical |
| **Target extraction** | 25 antigens detected: CD19, BCMA, CD22, CD20, CD30, HER2, GPC3, etc. |

### 3.2 `cart_trials` — 973 records

| Attribute | Value |
|---|---|
| **Source** | ClinicalTrials.gov REST API v2 |
| **Ingest time** | ~3 min |
| **Fields** | NCT ID, title, text_summary, phase, status, sponsor, target_antigen, car_generation, costimulatory, disease, enrollment, start_year, outcome_summary |
| **Phase distribution** | Early Phase 1 through Phase 3 |
| **Status** | Recruiting, completed, terminated, withdrawn, active |
| **Antigen extraction** | Automated from trial title and description |

### 3.3 `cart_constructs` — 6 records

| Attribute | Value |
|---|---|
| **Source** | FDA-approved CAR-T product labels |
| **Products** | Kymriah (tisagenlecleucel), Yescarta (axicabtagene ciloleucel), Tecartus (brexucabtagene autoleucel), Breyanzi (lisocabtagene maraleucel), Abecma (idecabtagene vicleucel), Carvykti (ciltacabtagene autoleucel) |
| **Fields** | name, text_summary, target_antigen, scfv_origin, costimulatory_domain, signaling_domain, generation, hinge_tm, vector_type, fda_status, known_toxicities |
| **Construct IDs** | `fda-tisagenlecleucel`, `fda-axicabtagene-ciloleucel`, `fda-brexucabtagene-autoleucel`, `fda-lisocabtagene-maraleucel`, `fda-idecabtagene-vicleucel`, `fda-ciltacabtagene-autoleucel` |

### 3.4 `cart_assays` — 45 records

| Attribute | Value |
|---|---|
| **Source** | Curated from landmark publications |
| **Papers** | ELIANA (NEJM 2018), ZUMA-1 (NEJM 2017), ZUMA-2 (NEJM 2020), TRANSFORM (Lancet 2022), KarMMa (NEJM 2021), CARTITUDE-1 (Lancet 2021), plus preclinical studies |
| **Assay types** | Cytotoxicity (12), in vivo/clinical (9), persistence (5), flow cytometry (5), exhaustion (3), cytokine (3), proliferation (3) |
| **Coverage** | All 6 FDA products + CD22, dual-target, GPC3, HER2, Mesothelin |
| **Resistance data** | CD19 loss/mutation, trogocytosis, lineage switch, BCMA biallelic loss, sBCMA decoy |
| **Linked constructs** | Records reference FDA construct IDs where applicable |

### 3.5 `cart_manufacturing` — 30 records

| Attribute | Value |
|---|---|
| **Source** | Curated from published manufacturing data and FDA guidance |
| **Process steps** | Transduction (6), expansion (7), harvest (2), formulation (3), cryopreservation (2), release testing (6), emerging platforms (4) |
| **Vector types** | Lentiviral, gamma-retroviral, transposon (Sleeping Beauty/piggyBac), mRNA |
| **Coverage** | VCN, titer, transduction efficiency, IL-2 vs IL-7/IL-15 expansion, rapid 6-day (Kite), defined CD4:CD8 (Breyanzi), POC manufacturing, allogeneic/off-the-shelf, cost analysis, GMP facility requirements |
| **Key parameters** | Functional titer, VCN, transduction efficiency, fold expansion, viability, sterility, CAR expression, RCL, identity |

### 3.6 Index Configuration (all collections)

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 1024 (literature), 256 (trials), 128 (constructs, assays, manufacturing) |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Knowledge Graph

### 4.1 Target Antigens (25 entries)

Each entry includes: protein name, UniProt ID, expression pattern, associated diseases, FDA-approved products, key clinical trials, known resistance mechanisms, toxicity profile, and normal tissue expression.

| Target | Diseases | Approved Products |
|---|---|---|
| **CD19** | B-ALL, DLBCL, FL, MCL, CLL | Kymriah, Yescarta, Tecartus, Breyanzi |
| **BCMA** | Multiple Myeloma | Abecma, Carvykti |
| CD22 | B-ALL (CD19-neg relapse) | — |
| CD20 | NHL, CLL | — |
| CD30 | Hodgkin lymphoma | — |
| CD33 | AML | — |
| CD38 | Multiple Myeloma | — |
| CD123 | AML, BPDCN | — |
| GD2 | Neuroblastoma | — |
| HER2 | Breast, gastric (solid tumor) | — |
| GPC3 | Hepatocellular carcinoma | — |
| EGFR | Glioblastoma, NSCLC | — |
| Claudin18.2 | Gastric, pancreatic | — |
| Mesothelin | Mesothelioma, ovarian, pancreatic | — |
| + 11 more | Various | — |

### 4.2 Toxicity Profiles (8 entries)

| Toxicity | Mechanism | Management |
|---|---|---|
| **CRS** | IL-6/IFN-γ cytokine storm | Tocilizumab, corticosteroids |
| **ICANS** | CNS endothelial activation | Dexamethasone, supportive |
| **B-cell aplasia** | On-target CD19/CD22 depletion | IVIG replacement |
| **HLH/MAS** | Macrophage hyperactivation | Anakinra, etoposide |
| **Cytopenias** | Marrow suppression | G-CSF, transfusions |
| **TLS** | Rapid tumor lysis | Rasburicase, hydration |
| **GvHD** | Allogeneic only, donor T-cells | Steroids, ruxolitinib |
| **On-target/off-tumor** | Normal tissue expression | Affinity tuning, safety switches |

### 4.3 Manufacturing Processes (10 entries)

Lentiviral transduction, retroviral transduction, T-cell expansion, leukapheresis, cryopreservation, release testing, vector production, quality control, formulation, potency testing.

### 4.4 Entity Aliases (39+ entries)

For Comparative Analysis Mode, the knowledge graph includes entity aliases that resolve product names, costimulatory domains, vector types, biomarkers, and regulatory terms to canonical entities with associated target antigens.

| Alias Category | Count | Examples |
|---|---|---|
| FDA Products | 12 | Kymriah → CD19, Carvykti → BCMA, tisagenlecleucel → CD19, etc. |
| Costimulatory Domains | 4 | 4-1BB (CD137), CD28, 4-1BB/CD28, OX40 |
| Vector Types | 2 | Lentiviral, Retroviral |
| Biomarker Terms | 8 | Ferritin, CRP, IL-6, sIL-2R, LDH, etc. |
| Regulatory Terms | 6 | BLA, RMAT, accelerated approval, breakthrough therapy, etc. |
| Safety Terms | 7 | REMS, FAERS, black box warning, post-marketing, etc. |

**Resolution priority:** CART_TARGETS (25) → ENTITY_ALIASES (39+) → CART_TOXICITIES (8) → CART_MANUFACTURING (10)

### 4.5 API Functions

```python
get_target_context("CD19")              # Returns full CD19 knowledge block
get_toxicity_context("CRS")             # Returns CRS management details
get_manufacturing_context(...)           # Returns manufacturing process details
get_all_context_for_query(text)          # Auto-detects entities and returns all relevant context
get_knowledge_stats()                    # Returns counts: {target_antigens: 25, ...}
resolve_comparison_entity("Kymriah")     # → {"type": "product", "canonical": "Kymriah (tisagenlecleucel)", "target": "CD19"}
get_comparison_context(entity_a, entity_b)  # Side-by-side knowledge graph context for dual entities
```

---

## 5. Query Expansion

Twelve expansion map categories (6 original + 6 added for expanded collections):

| Category | Keywords | Expanded Terms | Examples |
|---|---|---|---|
| Target Antigen | 26 | 196 | CD19 → [CD19, B-ALL, DLBCL, tisagenlecleucel, axicabtagene, ...] |
| Disease | 16 | 143 | multiple myeloma → [MM, plasma cell neoplasm, RRMM, ...] |
| Toxicity | 14 | 136 | CRS → [cytokine release syndrome, cytokine storm, tocilizumab, IL-6, ...] |
| Manufacturing | 16 | 181 | transduction → [lentiviral, retroviral, VCN, MOI, viral vector, ...] |
| Mechanism | 19 | 224 | resistance → [antigen loss, lineage switch, trogocytosis, exhaustion, ...] |
| Construct | 20 | 206 | bispecific → [dual-targeting, tandem, bicistronic, CD19/CD22, ...] |
| Safety | 15 | 135 | REMS → [risk evaluation, mitigation strategy, FAERS, adverse event, ...] |
| Biomarker | 14 | 125 | CRS prediction → [ferritin, CRP, IL-6, sIL-2R, predictive biomarker, ...] |
| Regulatory | 12 | 108 | BLA → [biologics license application, accelerated approval, RMAT, ...] |
| Sequence | 10 | 95 | scFv → [single-chain variable fragment, VH, VL, CDR, binding affinity, ...] |
| Real-World | 12 | 110 | registry → [real-world evidence, CIBMTR, post-marketing, outcomes, ...] |
| Immunogenicity | 10 | 92 | ADA → [anti-drug antibody, immunogenicity, neutralizing antibody, ...] |
| **Total** | **169** | **1,496** | |

The `expand_query()` function detects keywords in the user's query and returns relevant expansion terms, which are used to run additional filtered searches across collections.

---

## 6. Multi-Collection RAG Engine

### 6.1 Search Flow (measured on DGX Spark)

```
User Query: "Why do CD19 CAR-T therapies fail in relapsed B-ALL?"
    │
    ├── 1. Embed query with BGE asymmetric prefix               [< 5 ms]
    │      "Represent this sentence for searching relevant passages: ..."
    │
    ├── 2. Parallel search across 11 collections (top-5 each)   [12-16 ms]
    │   ├── cart_literature:     CD19 CAR-T failure papers       (score: 0.82-0.90)
    │   ├── cart_trials:         Terminated CD19 B-ALL trials    (score: 0.74-0.85)
    │   ├── cart_constructs:     CD19 CAR designs                (score: 0.78-0.86)
    │   ├── cart_assays:         Resistance/failure assay data   (score: 0.76-0.85)
    │   └── cart_manufacturing:  Production failure modes        (score: 0.71-0.79)
    │
    ├── 3. Query expansion: "CD19" → [CD19, B-ALL, DLBCL,       [< 1 ms]
    │      tisagenlecleucel, axicabtagene, ...]
    │
    ├── 4. Expanded filtered search (top-3 per term, top-5       [8-12 ms]
    │      expansion terms, collections with target_antigen)
    │
    ├── 5. Merge + deduplicate + rank (cap at 30 results)        [< 1 ms]
    │      Weights: lit 0.30, trial 0.25, construct 0.20,
    │               assay 0.15, manufacturing 0.10
    │
    ├── 6. Knowledge graph augmentation:                         [< 1 ms]
    │      CD19 → known_resistance: [CD19 loss, lineage switch, trogocytosis]
    │      CD19 → toxicity_profile: {CRS: 30-90%, ICANS: 20-65%}
    │      CD19 → approved_products: [Kymriah, Yescarta, Tecartus, Breyanzi]
    │
    ├── 7. Build prompt: evidence grouped by collection +         [< 1 ms]
    │      knowledge context + question + citation instructions
    │
    └── 8. Stream Claude Sonnet 4.6 response                    [~22-24 sec]
           Grounded cross-functional answer with clickable
           PubMed and ClinicalTrials.gov citation links
```

**Total: ~24 sec** (dominated by LLM generation; retrieval is ~25 ms)

### 6.2 Collection Weights

| Collection | Weight | Rationale |
|---|---|---|
| cart_literature | 0.30 | Published evidence is the primary source of truth |
| cart_trials | 0.25 | Clinical outcomes provide direct translational answers |
| cart_constructs | 0.20 | Design data explains mechanisms and structure-function |
| cart_assays | 0.15 | Lab data supports mechanistic claims with quantitative evidence |
| cart_manufacturing | 0.10 | Manufacturing links to clinical outcomes and feasibility |

### 6.3 System Prompt

The agent uses a specialized system prompt instructing Claude to:
1. **Cite evidence with clickable links** — Literature citations link to PubMed (`[Literature:PMID](https://pubmed.ncbi.nlm.nih.gov/PMID/)`), trial citations link to ClinicalTrials.gov (`[Trial:NCT...](https://clinicaltrials.gov/study/NCT...)`)
2. **Think cross-functionally** — connect insights across development stages
3. **Highlight failure modes** and resistance mechanisms when relevant
4. **Be specific** — cite trial names (ELIANA, ZUMA-1), products (Kymriah, Yescarta), and quantitative data
5. **Acknowledge uncertainty** — distinguish established facts from emerging data
6. **Suggest optimization strategies** based on historical data

### 6.4 Embedding Strategy

BGE-small-en-v1.5 uses asymmetric encoding — queries and documents are embedded differently:

| Mode | Prefix | Usage |
|---|---|---|
| **Query** | `"Represent this sentence for searching relevant passages: "` | User questions via `_embed_query()` |
| **Document** | None (raw text) | Ingested records via `to_embedding_text()` |

This asymmetric approach improves retrieval relevance by 5-15% compared to symmetric encoding.

### 6.5 Comparative Analysis Mode

Comparative queries are **auto-detected** and produce structured side-by-side analysis with markdown tables, advantages/limitations, and clinical context.

#### Detection and Parsing

```
User Query: "Compare 4-1BB vs CD28 costimulatory domains"
    │
    ├── 1. _is_comparative() detects COMPARE/VS/VERSUS/COMPARING       [< 1 ms]
    │
    ├── 2. _parse_comparison_entities() extracts raw entities            [< 1 ms]
    │      Pattern 1: "X vs/versus Y" (greedy group 2)
    │      Pattern 2: "Compare X and/with Y"
    │      Post-processing: strip trailing context words
    │        (costimulatory domains, resistance mechanisms, etc.)
    │
    ├── 3. resolve_comparison_entity() for each raw entity               [< 1 ms]
    │      Priority: CART_TARGETS → ENTITY_ALIASES → CART_TOXICITIES → CART_MANUFACTURING
    │      "4-1BB" → {"type": "costimulatory", "canonical": "4-1BB (CD137)"}
    │      "CD28"  → {"type": "costimulatory", "canonical": "CD28"}
    │
    ├── 4. Dual retrieve() — one per entity                              [~365 ms]
    │      Entity A: retrieve(question, target_antigen=entity_a.target)
    │      Entity B: retrieve(question, target_antigen=entity_b.target)
    │      ~24 results per entity across 11 collections
    │
    ├── 5. get_comparison_context() — side-by-side knowledge graph       [< 1 ms]
    │      Calls get_target_context() / get_toxicity_context() /
    │      get_manufacturing_context() for each entity
    │
    ├── 6. _build_comparative_prompt() — structured prompt               [< 1 ms]
    │      Evidence grouped by entity (A section, B section)
    │      Knowledge context appended
    │      Instructions: comparison table + advantages + limitations
    │
    └── 7. Stream Claude Sonnet 4.6 (max_tokens=3000)                   [~28-30 sec]
           Structured output: table, pros/cons, clinical context
```

**Total: ~30 sec** (365ms retrieval + ~30 sec LLM generation)

#### Supported Entity Types

| Entity Type | Examples | Resolution |
|---|---|---|
| Target Antigens | CD19, BCMA, CD22, CD20 | Direct match in CART_TARGETS (25 entries) |
| FDA Products | Kymriah, Yescarta, Carvykti, Abecma | ENTITY_ALIASES → canonical name + target |
| Costimulatory Domains | 4-1BB, CD28 | ENTITY_ALIASES → type: costimulatory |
| Toxicity Profiles | CRS, ICANS, B-cell aplasia | CART_TOXICITIES → type: toxicity |
| Manufacturing Processes | Lentiviral, transduction | CART_MANUFACTURING → type: manufacturing |

#### Example Comparative Queries

```
"Compare CD19 vs BCMA"                              → Target vs target (with target_antigen filtering)
"Compare 4-1BB vs CD28 costimulatory domains"        → Costimulatory domain comparison
"Kymriah versus Carvykti"                            → Product vs product (resolves to CD19 vs BCMA)
"Compare CRS and ICANS toxicity"                     → Toxicity profile comparison
```

#### Fallback Behavior

If entity parsing fails (unrecognized entities, ambiguous input), the query gracefully falls back to the standard single-query `retrieve()` path with no user-visible error.

#### UI Evidence Panel

Comparative evidence is displayed in an entity-grouped collapsible panel:

- **Entity A header** (blue) — evidence cards for entity A results
- **"— VS —" divider** (green)
- **Entity B header** (purple) — evidence cards for entity B results

Each evidence card shows collection badge, ID, cosine score, clickable source link, and text snippet.

### 6.6 Evidence Panel and Clickable Citations

All query responses include a collapsible evidence panel with collection-badged cards:

| Badge Color | Collection | Link Format |
|---|---|---|
| Blue | Literature | [PubMed](https://pubmed.ncbi.nlm.nih.gov/PMID/) |
| Purple | Trial | [ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT...) |
| Green | Construct | — |
| Yellow | Assay | — |
| Orange | Manufacturing | — |

Each card displays: collection badge, record ID, cosine similarity score, clickable source link (where available), and a 200-character text snippet.

---

## 7. Data Sources and Ingest Pipelines

### 7.1 Actual Ingest Results

| Source | API | Records | Time | Collection |
|---|---|---|---|---|
| PubMed | NCBI E-utilities (esearch + efetch) | 5,047 | ~15 min | cart_literature |
| ClinicalTrials.gov | REST API v2 | 973 | ~3 min | cart_trials |
| FDA Product Labels | Manual curation (seed script) | 6 | ~5 sec | cart_constructs |
| Published Assays | Curated JSON (seed script) | 45 | ~30 sec | cart_assays |
| Published Manufacturing | Curated JSON (seed script) | 30 | ~30 sec | cart_manufacturing |
| Safety | Curated pharmacovigilance data | 40 | ~30 sec | cart_safety |
| Biomarkers | Curated CRS/exhaustion markers | 43 | ~30 sec | cart_biomarkers |
| Regulatory | FDA approval timelines | 25 | ~30 sec | cart_regulatory |
| Sequences | scFv/molecular binding data | 27 | ~30 sec | cart_sequences |
| Real-World Evidence | Registry outcomes | 30 | ~30 sec | cart_realworld |
| **Total** | | **6,266** | **~22 min** | |

### 7.2 Ingest Pipeline Architecture

All 5 ingest pipelines inherit from `BaseIngestPipeline` with a standard flow:

```
fetch(**kwargs)           # Retrieve raw data (API call, file read)
    │
    ▼
parse(raw_data)           # Validate into Pydantic models
    │
    ▼
embed_and_store(records)  # Embed text → insert into Milvus
    │
    ├── record.to_embedding_text()   # Generate embedding input
    ├── embedder.encode(texts)       # BGE-small batch encoding
    ├── Enum → str conversion        # ProcessStep, AssayType, etc.
    ├── UTF-8 byte truncation        # Milvus VARCHAR byte limits
    └── manager.insert_batch()       # Batch insert into collection
```

### 7.3 Assay Seed Data Coverage

| Category | Records | Key Data |
|---|---|---|
| Cytotoxicity | 12 | All 6 FDA products, CD22, dual-target, GPC3, HER2, Mesothelin, sBCMA decoy |
| In vivo / Clinical | 9 | ELIANA, ZUMA-1, ZUMA-2, TRANSFORM, KarMMa, CARTITUDE-1, mouse models |
| Persistence | 5 | 4-1BB vs CD28 persistence kinetics for all products |
| Flow cytometry | 5 | CD19 loss, trogocytosis, lineage switch, BCMA loss, product characterization |
| Exhaustion | 3 | CD28 vs 4-1BB exhaustion, ZUMA-1/KarMMa correlative data |
| Cytokine | 3 | IFN-gamma profiles: tisagenlecleucel, axicabtagene, lisocabtagene |
| Proliferation | 3 | IL-7/IL-15 vs IL-2, tisagenlecleucel 42x, axicabtagene 25x (6-day rapid) |

### 7.4 Manufacturing Seed Data Coverage

| Category | Records | Key Data |
|---|---|---|
| Transduction | 6 | Lentiviral titer/VCN/efficiency, retroviral, transposon, mRNA |
| Expansion | 7 | IL-2, IL-7/IL-15, rapid 6-day (Kite), defined CD4:CD8 (BMS), POC, T-cell selection, failure modes |
| Harvest / Formulation | 4 | Leukapheresis yield, viability, dosing, vein-to-vein time |
| Cryopreservation | 2 | Controlled-rate freezing, shipping logistics (LN2 dry shippers) |
| Release testing | 6 | Sterility (USP <71>), CAR expression, RCL/RCR, identity/COI, potency, residual beads |
| Emerging / Cost | 5 | Allogeneic off-the-shelf, cost breakdown ($150-300K), GMP facility, lymphodepletion, capacity |

---

## 8. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 8.1 Ingest Performance

| Operation | Time | Records | Rate |
|---|---|---|---|
| PubMed fetch + parse + embed + store | ~15 min | 5,047 | ~5.6 rec/sec |
| ClinicalTrials.gov fetch + embed + store | ~3 min | 973 | ~5.4 rec/sec |
| Assay seed embed + store | ~30 sec | 45 | ~1.5 rec/sec |
| Manufacturing seed embed + store | ~30 sec | 30 | ~1.0 rec/sec |
| FDA construct seed (6 products) | ~5 sec | 6 | ~1.2 rec/sec |

Note: Ingest rate is dominated by BGE-small embedding time (~180ms per text on CPU). GPU acceleration would increase throughput 10-50x.

### 8.2 Search Performance

| Operation | Latency | Notes |
|---|---|---|
| Single collection search (top-5) | 3-5 ms | Milvus IVF_FLAT with cached index |
| 11-collection parallel search (top-5 each) | 12-16 ms | Sequential per-collection, 55 total results |
| Query expansion + filtered search | 8-12 ms | Up to 5 expanded terms × applicable collections |
| Knowledge graph augmentation | < 1 ms | In-memory dictionary lookup |
| Full retrieve() pipeline | 20-30 ms | Embed + search + expand + merge + knowledge |
| Comparative entity resolution | < 1 ms | CART_TARGETS → ENTITY_ALIASES → TOXICITIES → MFG |
| Comparative dual retrieval (2 × 11 collections) | ~365 ms | Two retrieve() calls, ~46 total results (24 + 22) |

### 8.3 RAG Query Performance

| Operation | Latency | Notes |
|---|---|---|
| Full query() (retrieve + Claude generate) | ~24 sec | Dominated by LLM generation |
| Comparative query (dual retrieve + Claude) | ~30 sec | 365ms retrieval + structured comparison prompt |
| Streaming query_stream() (time to first token) | ~3 sec | Evidence returned immediately |
| Response length (standard) | 800-2000 chars | Grounded answer with citations |
| Response length (comparative) | 1500-3000 chars | Structured tables + pros/cons + clinical context |
| Token count (standard) | 400-800 tokens | Claude Sonnet 4.6 output |
| Token count (comparative) | 800-1500 tokens | Claude Sonnet 4.6 output (max_tokens=3000) |

### 8.4 Relevance Quality

| Query | Top Hit Score | Collection | Expected? |
|---|---|---|---|
| "CD19 CAR-T cytotoxicity" | 0.791 | cart_assays | Yes (CD22/CD19 assay data) |
| "BCMA resistance mechanisms" | 0.809 | cart_assays | Yes (sBCMA decoy mechanism) |
| "lentiviral VCN transduction" | 0.870 | cart_manufacturing | Yes (VCN quality attribute) |
| "CAR-T shipping cryopreservation" | 0.779 | cart_manufacturing | Yes (cryo shipping logistics) |
| "CD19 CAR-T failure B-ALL" | 0.82-0.90 | cart_literature | Yes (PubMed abstracts) |

---

## 9. Infrastructure

### 9.1 Technology Stack

| Component | Technology | Version/Detail |
|---|---|---|
| Vector database | Milvus | 2.4, localhost:19530 |
| Embedding model | BGE-small-en-v1.5 | 384-dim, BAAI, ~33M params |
| LLM | Claude Sonnet 4.6 | Anthropic API, claude-sonnet-4-20250514 |
| UI framework | Streamlit | Port 8521, NVIDIA black/green theme |
| Data models | Pydantic | BaseModel + Field validation |
| Configuration | Pydantic BaseSettings | Environment variable support |
| Hardware target | NVIDIA DGX Spark | GB10 GPU, 128GB unified, $3,999 |

### 9.2 Service Ports

| Port | Service |
|---|---|
| 8521 | CAR-T Intelligence Agent Streamlit UI |
| 19530 | Milvus vector database (shared with main pipeline) |

### 9.3 Dependencies on HCLS AI Factory

| Dependency | Usage |
|---|---|
| Milvus 2.4 instance | Shared vector database — CAR-T adds 10 owned collections alongside existing `genomic_evidence` (read-only) |
| `ANTHROPIC_API_KEY` | Loaded from `rag-chat-pipeline/.env` if not set in environment |
| BGE-small-en-v1.5 | Same embedding model as main RAG pipeline |

---

## 10. Demo Scenarios

### 10.1 Validated Demo Queries

These queries have been tested end-to-end through the RAG pipeline:

**1. "Why do CD19 CAR-T therapies fail in relapsed B-ALL patients?"**
- Searches: literature (resistance papers), trials (terminated CD19 trials), assays (CD19 loss, trogocytosis, lineage switch data), constructs (CD19 product designs)
- Knowledge graph: CD19 → known_resistance, toxicity_profile, approved_products
- Expected answer: Covers antigen loss (28% of relapses), lineage switch (10%, KMT2A-associated), trogocytosis, T-cell exhaustion

**2. "Compare 4-1BB vs CD28 costimulatory domains for DLBCL"**
- Searches: literature (head-to-head reviews), trials (ZUMA-1 vs TRANSCEND/TRANSFORM), assays (exhaustion markers, persistence data), constructs (Yescarta vs Breyanzi)
- Expected answer: CD28 = faster kinetics, higher peak, more exhaustion; 4-1BB = sustained persistence, lower toxicity

**3. "What manufacturing parameters predict clinical response?"**
- Searches: literature (correlative studies), manufacturing (VCN, expansion, phenotype), assays (product characterization), trials (responder analyses)
- Expected answer: T-cell fitness, Tcm frequency (>40% threshold), CD4:CD8 ratio, VCN, manufacturing time

**4. "BCMA CAR-T resistance mechanisms in multiple myeloma"**
- Searches: literature (resistance reviews), assays (biallelic BCMA loss, sBCMA decoy), constructs (ide-cel vs cilta-cel), trials (KarMMa/CARTITUDE relapse data)
- Expected answer: BCMA downregulation, biallelic loss (29% of relapses), soluble BCMA decoy, antigen density heterogeneity

**5. "How does T-cell exhaustion affect CAR-T persistence?"**
- Searches: literature (exhaustion biology), assays (PD-1/LAG-3/TIM-3 data, CD28 vs 4-1BB exhaustion), manufacturing (IL-7/IL-15 vs IL-2)
- Expected answer: Exhaustion markers predict outcomes, CD28 drives faster exhaustion, 4-1BB preserves Tcm, IL-7/IL-15 expansion reduces exhaustion

---

## 11. File Structure (Actual)

```
cart_intelligence_agent/
├── Docs/
│   └── CART_Intelligence_Agent_Design.md    # This document
├── src/
│   ├── __init__.py
│   ├── models.py                            # Pydantic data models + ComparativeResult (299 lines)
│   ├── collections.py                       # Milvus collection schemas + manager (842 lines)
│   ├── knowledge.py                         # Knowledge graph + entity aliases + comparison (1,030 lines)
│   ├── query_expansion.py                   # 12 expansion maps, 169→1496 terms (955 lines)
│   ├── rag_engine.py                        # Multi-collection RAG + comparative analysis (558 lines)
│   ├── agent.py                             # CAR-T Intelligence Agent (262 lines)
│   ├── ingest/
│   │   ├── __init__.py
│   │   ├── base.py                          # Base ingest pipeline (184 lines)
│   │   ├── literature_parser.py             # PubMed E-utilities ingest (350 lines)
│   │   ├── clinical_trials_parser.py        # ClinicalTrials.gov API v2 ingest (403 lines)
│   │   ├── construct_parser.py              # CAR construct parser + FDA seed (292 lines)
│   │   ├── assay_parser.py                  # Assay data parser (163 lines)
│   │   └── manufacturing_parser.py          # Manufacturing/CMC parser (111 lines)
│   └── utils/
│       ├── __init__.py
│       └── pubmed_client.py                 # NCBI E-utilities HTTP client (390 lines)
├── app/
│   └── cart_ui.py                           # Streamlit chat + comparative UI (484 lines)
├── config/
│   └── settings.py                          # Pydantic BaseSettings (69 lines)
├── data/
│   ├── reference/
│   │   ├── assay_seed_data.json             # 45 curated assay records
│   │   └── manufacturing_seed_data.json     # 30 curated manufacturing records
│   └── cache/                               # Embedding cache
├── scripts/
│   ├── setup_collections.py                 # Create collections + seed FDA constructs (86 lines)
│   ├── ingest_pubmed.py                     # CLI: PubMed ingest (185 lines)
│   ├── ingest_clinical_trials.py            # CLI: ClinicalTrials.gov ingest (182 lines)
│   ├── seed_assays.py                       # CLI: Seed assay data (95 lines)
│   ├── seed_manufacturing.py                # CLI: Seed manufacturing data (95 lines)
│   ├── validate_e2e.py                      # End-to-end 5-test validation (191 lines)
│   ├── test_rag_pipeline.py                 # Full RAG + LLM integration test (223 lines)
│   └── seed_knowledge.py                    # Knowledge graph export (110 lines)
├── requirements.txt
├── LICENSE                                  # Apache 2.0
└── README.md
```

**55 Python files | ~16,748 lines of code | Apache 2.0**

---

## 12. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Scaffold** | Complete | Data models, collection schemas, knowledge graph, query expansion, RAG engine, agent, Streamlit UI, ingest pipeline stubs |
| **Phase 2: Data Ingest** | Complete | PubMed (5,047), ClinicalTrials.gov (973), FDA constructs (6), assay seed (45), manufacturing seed (30), safety (40), biomarkers (43), regulatory (25), sequences (27), real-world (30) |
| **Phase 3: RAG Integration** | Complete | Multi-collection search, knowledge augmentation, Claude Sonnet 4.6 streaming, 5 demo queries validated |
| **Phase 4: UI + Demo** | Complete | Streamlit UI on port 8521, NVIDIA theme, sidebar filters, demo query buttons, streaming responses |
| **Phase 5: UI + Analysis** | Complete | Clickable PubMed/ClinicalTrials.gov citation links, collapsible evidence panel with collection badges, **Comparative Analysis Mode** with auto-detection, entity resolution (39+ aliases), dual retrieval (~365ms), entity-grouped evidence panel, and structured markdown comparison tables |

### Remaining Work

| Item | Priority | Effort |
|---|---|---|
| Export to report (PDF/markdown summary generation) | Medium | 2-3 hours |
| Agent reasoning loop testing (`agent.py` plan→search→synthesize) | Medium | 1-2 hours |
| Genomic evidence bridge (connect to existing `genomic_evidence` collection) | Low | 2-3 hours |
| Nextflow orchestrator integration | Low | 1-2 hours |
| Additional construct data (published designs beyond FDA-approved) | Low | 2-3 hours |
| Landing page integration (health check endpoint) | Low | 1 hour |

---

## 13. Relationship to HCLS AI Factory

This agent demonstrates the **generalizability** of the HCLS AI Factory architecture. The same infrastructure that supports the VCP/Frontotemporal Dementia pipeline extends to CAR-T cell therapy intelligence with:

- **Same Milvus instance** — 10 new owned collections alongside existing `genomic_evidence` (3.56M vectors, read-only)
- **Same embedding model** — BGE-small-en-v1.5 (384-dim)
- **Same LLM** — Claude via Anthropic API
- **Same hardware** — NVIDIA DGX Spark ($3,999)
- **Same patterns** — Pydantic models, BaseIngestPipeline, knowledge graph, query expansion

The key architectural insight: **the platform is not disease-specific**. By changing the knowledge graph, query expansion maps, and collection schemas, the same RAG architecture serves any therapeutic area. The CAR-T agent proves this with a completely different domain (cell therapy manufacturing vs. small molecule drug discovery) running on the same infrastructure.

---

## 14. Credits

- **Adam Jones** — HCLS AI Factory platform, 14+ years genomic research
- **Apache 2.0 License**
