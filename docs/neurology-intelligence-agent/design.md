# Neurology Intelligence Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## 1. Executive Summary

The Neurology Intelligence Agent extends the HCLS AI Factory platform to deliver RAG-powered clinical decision support across the full spectrum of neurological disease. It unifies fragmented neurological evidence — spanning cerebrovascular, neurodegenerative, epilepsy, movement disorders, multiple sclerosis, headache, neuromuscular, and neuro-oncology domains — into a single intelligence platform. Clinicians receive guideline-grounded, evidence-cited recommendations with validated clinical scale calculations in under five seconds.

The system implements **10 validated clinical scale calculators** (NIHSS, GCS, MoCA, EDSS, Hoehn-Yahr, MDS-UPDRS Part III, HIT-6, ALSFRS-R, ASPECTS, mRS), provides **8 evidence-based clinical workflows**, and searches **14 Milvus vector collections** containing disease-specific literature, clinical trials, imaging protocols, treatment guidelines, and genomic correlations.

The platform enables time-critical queries like *"Acute left MCA stroke, NIHSS 18, last known well 3 hours ago — treatment options?"* that simultaneously search stroke protocols, imaging guidelines, treatment evidence, and clinical trials — returning grounded recommendations with guideline citations and scale calculations.

### Key Results

| Metric | Value |
|---|---|
| Milvus collections | **14** domain-specific collections (13 owned + 1 read-only) |
| Clinical scale calculators | **10** validated instruments |
| Clinical workflows | **8** + 1 general neurological query |
| Disease domains | **10** (stroke, dementia, epilepsy, brain tumor, MS, Parkinson's, headache, neuromuscular, neurocritical, neuropathy) |
| Query expansion aliases | **251+** across 16 synonym maps |
| Estimated total records | **~855,000** across all collections |
| Test suite | **209** automated tests across 12 modules |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Neurology Agent Role |
|---|---|
| **DataStore** | Raw files: PubMed XML, ClinicalTrials.gov JSON, guideline documents, imaging protocol specs |
| **DataEngine** | Ingest pipelines for literature, trials, guidelines, imaging, drug, genomic, and pathway data |
| **DataBase** | 14 Milvus collections (13 owned + 1 read-only) + knowledge base (disease taxonomy, drug catalog, guideline recommendations) |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + 10 clinical scale calculators + query expansion (251+ aliases) |
| **AgentEngine** | NeurologyAgent orchestrator + Streamlit UI + FastAPI REST |

### 2.2 System Diagram

```
+================================================================+
|                     PRESENTATION LAYER                          |
|  +---------------------+   +-----------------------------+     |
|  | Streamlit Chat UI   |   | FastAPI REST API             |     |
|  | Port 8534           |   | Port 8528                    |     |
|  | Interactive Q&A     |   | Versioned endpoints (v1)     |     |
|  | Scale calculators   |   | CORS, auth, rate limiting    |     |
|  +----------+----------+   +-------------+---------------+     |
+================================================================+
               |                            |
+================================================================+
|                    INTELLIGENCE LAYER                           |
|  +----------------+  +----------------+  +------------------+  |
|  | Query Expander |  | Workflow Engine |  | Clinical Scale   |  |
|  | 251+ aliases   |  | 8+1 workflows  |  | Calculators      |  |
|  | 16 synonym     |  | domain-specific|  | 10 validated     |  |
|  | maps           |  | weight boost   |  | instruments      |  |
|  +-------+--------+  +-------+--------+  +--------+---------+  |
|          |                   |                     |            |
|  +-------v-------------------v---------------------v---------+  |
|  |              Neurology RAG Engine                         |  |
|  | ThreadPoolExecutor parallel search across 14 collections  |  |
|  | Workflow-specific collection weight boosting              |  |
|  | Citation scoring (high/medium/low)                        |  |
|  | Conversation memory (24h TTL)                             |  |
|  +---------------------------+-------------------------------+  |
|                              |                                  |
|  +---------------------------v-------------------------------+  |
|  |              Claude Sonnet 4.6 (Anthropic)                |  |
|  | Evidence synthesis with clinical system prompt            |  |
|  | Guideline-grounded recommendations (AAN/AHA/ILAE/MDS)    |  |
|  +-----------------------------------------------------------+  |
+================================================================+
               |
+================================================================+
|                        DATA LAYER                              |
|  +-----------------------------------------------------------+  |
|  |              Milvus 2.4 Vector Database                   |  |
|  |  14 collections | BGE-small 384-dim | IVF_FLAT/COSINE    |  |
|  |  855K estimated records | etcd + MinIO backend            |  |
|  +-----------------------------------------------------------+  |
+================================================================+
```

---

## 3. Data Collections — Actual State

### 3.1 Collection Catalog

| # | Collection | Est. Records | Weight | Primary Use |
|---|---|---|---|---|
| 1 | `neuro_literature` | 500,000 | 0.12 | PubMed neurological literature |
| 2 | `neuro_trials` | 50,000 | 0.09 | ClinicalTrials.gov neurology trials |
| 3 | `neuro_guidelines` | 2,000 | 0.10 | AAN/AHA/ILAE/MDS/IHS clinical guidelines |
| 4 | `neuro_imaging` | 5,000 | 0.08 | MRI/CT/PET protocols and findings |
| 5 | `neuro_drugs` | 3,000 | 0.08 | Neurological pharmacotherapy |
| 6 | `neuro_genomics` | 10,000 | 0.07 | Neurogenetic variants and associations |
| 7 | `neuro_pathways` | 2,000 | 0.06 | Neural pathways and circuits |
| 8 | `neuro_scales` | 500 | 0.06 | Clinical scale validation data |
| 9 | `neuro_surgery` | 3,000 | 0.05 | Neurosurgical procedures and outcomes |
| 10 | `neuro_electrophysiology` | 5,000 | 0.06 | EEG, EMG, NCS patterns and interpretation |
| 11 | `neuro_biomarkers` | 2,000 | 0.05 | CSF, serum, and imaging biomarkers |
| 12 | `neuro_rehabilitation` | 3,000 | 0.04 | Neurological rehabilitation protocols |
| 13 | `neuro_case_reports` | 5,000 | 0.05 | Published neurological case reports |
| 14 | `genomic_evidence` | ~265,000 | 0.03 | Shared genomic variant context |

### 3.2 Index Configuration

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 1024 (literature), 256 (trials), 128 (others) |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Clinical Scale Calculators

### 4.1 NIHSS (National Institutes of Health Stroke Scale)

| Domain | Items | Score Range |
|---|---|---|
| Level of consciousness | 3 items (LOC, LOC questions, LOC commands) | 0-7 |
| Gaze | 1 item | 0-2 |
| Visual fields | 1 item | 0-3 |
| Facial palsy | 1 item | 0-3 |
| Motor arm | 2 items (left, right) | 0-8 |
| Motor leg | 2 items (left, right) | 0-8 |
| Ataxia | 1 item | 0-2 |
| Sensory | 1 item | 0-2 |
| Language | 1 item | 0-3 |
| Dysarthria | 1 item | 0-2 |
| Extinction/inattention | 1 item | 0-2 |
| **Total** | **15 items** | **0-42** |
| **Output** | Score, severity (minor/moderate/moderate-severe/severe), tPA eligibility window |

### 4.2 GCS (Glasgow Coma Scale)

| Component | Range | Best Response |
|---|---|---|
| Eye opening | 1-4 | Spontaneous |
| Verbal response | 1-5 | Oriented |
| Motor response | 1-6 | Obeys commands |
| **Total** | **3-15** | Severity: mild (13-15), moderate (9-12), severe (3-8) |

### 4.3 MoCA (Montreal Cognitive Assessment)

| Domain | Max Points |
|---|---|
| Visuospatial/Executive | 5 |
| Naming | 3 |
| Memory | 5 (delayed recall) |
| Attention | 6 |
| Language | 3 |
| Abstraction | 2 |
| Orientation | 6 |
| **Total** | **30** (Normal >= 26, MCI 18-25, Dementia < 18) |

### 4.4 Additional Scales

| Scale | Domain | Range | Key Thresholds |
|---|---|---|---|
| **EDSS** | MS disability | 0-10 (0.5 steps) | 0 = normal, 6.0 = bilateral assistance, 10 = death |
| **Hoehn-Yahr** | Parkinson's staging | 1-5 | 1 = unilateral, 3 = bilateral with postural instability, 5 = wheelchair/bed |
| **MDS-UPDRS III** | PD motor exam | 0-132 | 18 items, each 0-4 severity scale |
| **HIT-6** | Headache impact | 36-78 | <= 49 little, 50-55 some, 56-59 substantial, >= 60 severe |
| **ALSFRS-R** | ALS function | 0-48 | 12 items: bulbar, fine motor, gross motor, respiratory |
| **ASPECTS** | Stroke CT | 0-10 | >= 7 favorable, < 7 large core |
| **mRS** | Disability outcome | 0-6 | 0 = no symptoms, 6 = dead |

---

## 5. Clinical Workflows

### 5.1 Workflow Catalog

| # | Workflow | Clinical Question | Key Scales | Weight-Boosted Collections |
|---|---|---|---|---|
| 1 | Acute Stroke | "tPA/thrombectomy eligibility for this stroke?" | NIHSS, ASPECTS, mRS | guidelines, imaging, drugs, literature |
| 2 | Dementia Evaluation | "ATN staging and anti-amyloid eligibility?" | MoCA, CDR | guidelines, biomarkers, genomics, drugs |
| 3 | Epilepsy Classification | "ILAE classification and drug-resistant assessment?" | Seizure frequency | guidelines, drugs, electrophysiology |
| 4 | Brain Tumor Grading | "WHO 2021 molecular classification?" | KPS | guidelines, genomics, imaging, surgery |
| 5 | MS Monitoring | "NEDA-3 status and DMT escalation?" | EDSS | guidelines, imaging, drugs, biomarkers |
| 6 | Parkinson's Assessment | "Motor severity and DBS candidacy?" | H&Y, UPDRS, MoCA | guidelines, drugs, surgery |
| 7 | Headache Classification | "ICHD-3 diagnosis and CGRP therapy guidance?" | HIT-6 | guidelines, drugs, imaging |
| 8 | Neuromuscular Evaluation | "ALS vs. neuropathy differential?" | ALSFRS-R | electrophysiology, genomics, drugs |

### 5.2 Workflow-Specific Weight Boosting

Each workflow dynamically adjusts collection weights. Example for Acute Stroke:

| Collection | Base Weight | Stroke Boost | Effective Weight |
|---|---|---|---|
| neuro_guidelines | 0.10 | 2.0x | 0.20 |
| neuro_imaging | 0.08 | 1.8x | 0.14 |
| neuro_drugs | 0.08 | 1.5x | 0.12 |
| neuro_literature | 0.12 | 1.0x | 0.12 |
| (others) | varies | 0.5-1.0x | reduced |

---

## 6. Multi-Collection RAG Engine

### 6.1 Search Flow

```
User Query: "Acute left MCA stroke, NIHSS 18, last known well 3h ago"
    │
    ├── 1. Workflow classification: ACUTE_STROKE                  [< 1 ms]
    │
    ├── 2. NIHSS calculation: Score 18 → Moderate-severe          [< 5 ms]
    │      tPA window: Within 4.5h → ELIGIBLE
    │      Thrombectomy: NIHSS >= 6, LVO suspected → EVALUATE
    │
    ├── 3. Embed query with BGE asymmetric prefix                 [< 5 ms]
    │
    ├── 4. Parallel search across 14 collections                  [12-18 ms]
    │      (with stroke workflow weight boosting)
    │   ├── neuro_guidelines:  AHA/ASA stroke guidelines (2x)    (score: 0.85-0.92)
    │   ├── neuro_imaging:     CTA/CTP/MRI stroke protocols (1.8x)(score: 0.80-0.88)
    │   ├── neuro_drugs:       tPA, TNK, antiplatelets (1.5x)    (score: 0.78-0.86)
    │   └── neuro_literature:  DAWN, DEFUSE-3 trials              (score: 0.75-0.84)
    │
    ├── 5. Query expansion: "MCA stroke NIHSS 18" →               [< 1 ms]
    │      [middle cerebral artery, large vessel occlusion,
    │       alteplase, tenecteplase, thrombectomy, DAWN, ...]
    │
    ├── 6. Knowledge base augmentation                            [< 1 ms]
    │
    └── 7. Stream Claude Sonnet 4.6 response                     [~22-26 sec]
           AHA/ASA guideline-grounded recommendation:
           tPA eligibility, thrombectomy criteria (DAWN/DEFUSE-3),
           imaging protocol, BP management
```

**Total: ~26 sec** (retrieval + scales: ~30 ms; LLM: ~25 sec)

### 6.2 Citation Scoring

| Level | Threshold | Display |
|---|---|---|
| High confidence | >= 0.75 | Full citation with source link |
| Medium confidence | >= 0.60 | Citation with caveat |
| Below threshold | < 0.40 | Filtered out |

---

## 7. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 7.1 Scale Calculator Performance

| Scale | Latency | Validated Against |
|---|---|---|
| NIHSS (15 items) | <5 ms | NIH Stroke Scale training |
| GCS (3 components) | <2 ms | Standard GCS protocol |
| MoCA (7 domains) | <5 ms | MoCA validation study |
| EDSS (functional systems) | <10 ms | Neurostatus-eEDSS |
| Hoehn-Yahr | <2 ms | Original Hoehn & Yahr 1967 |
| MDS-UPDRS III (18 items) | <10 ms | MDS-UPDRS validation |
| HIT-6 (6 items) | <2 ms | HIT-6 validation study |
| ALSFRS-R (12 items) | <5 ms | ALSFRS-R validation |
| ASPECTS (10 regions) | <5 ms | Original ASPECTS paper |
| mRS (single score) | <1 ms | Standard mRS protocol |
| **All 10 scales** | **<50 ms** | |

### 7.2 RAG Query Performance

| Operation | Latency |
|---|---|
| Full query (retrieve + Claude generate) | ~26 sec |
| Streaming query (time to first token) | ~3 sec |
| 14-collection parallel search | 12-18 ms |
| Query expansion | < 1 ms |
| Knowledge augmentation | < 1 ms |

---

## 8. Infrastructure

### 8.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4, localhost:19530 |
| Embeddings | BGE-small-en-v1.5 (BAAI) — 384-dim |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Web UI | Streamlit (port 8534, NVIDIA black/green theme) |
| REST API | FastAPI + Uvicorn (port 8528) |
| Configuration | Pydantic BaseSettings with `NEURO_` prefix |
| Testing | pytest (209 tests) |
| Hardware target | NVIDIA DGX Spark (GB10 GPU, 128GB unified, $4,699) |

### 8.2 Service Ports

| Port | Service |
|---|---|
| 8528 | FastAPI REST API |
| 8534 | Streamlit Chat UI |
| 19530 | Milvus vector database (shared) |

### 8.3 Dependencies on HCLS AI Factory

| Dependency | Usage |
|---|---|
| Milvus 2.4 instance | Shared vector database — adds 13 owned collections alongside existing `genomic_evidence` (read-only) |
| `ANTHROPIC_API_KEY` | Shared Anthropic API key |
| BGE-small-en-v1.5 | Same embedding model as main RAG pipeline |

---

## 9. Demo Scenarios

### 9.1 Validated Demo Queries

**1. "Acute left MCA stroke, NIHSS 18, last known well 3 hours ago — treatment?"**
- NIHSS: 18 → Moderate-severe stroke
- tPA: Within 4.5h window → Eligible (AHA/ASA Class I)
- Thrombectomy: NIHSS >= 6 + suspected LVO → Evaluate with CTA (DAWN/DEFUSE-3 criteria)

**2. "72-year-old with MoCA 22, amyloid PET positive — staging and treatment options?"**
- MoCA: 22 → Mild cognitive impairment
- ATN staging: A+ (amyloid PET) → Evaluate tau and neurodegeneration
- Anti-amyloid: Lecanemab eligibility criteria assessment

**3. "Drug-resistant epilepsy, failed 3 ASMs — surgical candidacy evaluation?"**
- ILAE definition: Failed 2+ appropriately chosen ASMs → Drug-resistant
- Surgical workup: Video-EEG monitoring, brain MRI (epilepsy protocol), neuropsych testing
- Evidence: Cochrane review on epilepsy surgery outcomes

**4. "WHO 2021 classification for IDH-mutant, 1p/19q codeleted brain tumor?"**
- Classification: Oligodendroglioma, IDH-mutant, 1p/19q-codeleted (WHO Grade 2 or 3)
- Treatment: PCV chemotherapy + radiation (RTOG 9402, EORTC 26951 evidence)

**5. "EDSS 4.0 MS patient on dimethyl fumarate with new enhancing lesions — escalation?"**
- EDSS: 4.0 → Ambulatory without aid, limited walking distance
- NEDA-3: Failed (new MRI activity)
- Escalation: Consider natalizumab, ocrelizumab, or ofatumumab (JCV status guides selection)

---

## 10. File Structure (Actual)

```
neurology_intelligence_agent/
├── src/
│   ├── agent.py                     # Agent orchestrator
│   ├── models.py                    # Enums and Pydantic models
│   ├── collections.py               # 14 collection schemas
│   ├── rag_engine.py                # Multi-collection RAG with weight boosting
│   ├── clinical_scales.py           # 10 validated scale calculators
│   ├── clinical_workflows.py        # 8 clinical workflows
│   ├── knowledge.py                 # Domain knowledge base
│   ├── query_expansion.py           # 251+ aliases, 16 synonym maps
│   ├── cross_modal.py               # Cross-agent triggers
│   ├── metrics.py                   # Prometheus metrics
│   └── export.py                    # Report formats
├── app/
│   └── neuro_ui.py                 # Streamlit chat interface
├── api/
│   └── main.py                     # FastAPI (15 clinical endpoints)
├── config/
│   └── settings.py                 # Pydantic BaseSettings (50+ params)
├── data/
│   ├── cache/                      # Conversation persistence (24h TTL)
│   └── reference/                  # Reference data files
├── scripts/
│   ├── setup_collections.py
│   ├── seed_knowledge.py
│   └── run_ingest.py
├── tests/                          # 209 tests across 12 modules
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```

---

## 11. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Architecture** | Complete | 14 collections, 10 scale calculators, 8 workflows, knowledge base, RAG engine |
| **Phase 2: Data** | Complete | 855K estimated records, 10 disease domains, guideline library (AAN/AHA/ILAE/MDS/IHS) |
| **Phase 3: RAG Integration** | Complete | Multi-collection parallel search with workflow-specific weight boosting, Claude streaming |
| **Phase 4: Testing** | Complete | 209 tests, all passing |
| **Phase 5: UI + Demo** | Complete | Streamlit UI on port 8534, NVIDIA theme, 5 demo scenarios validated |

### Remaining Work

| Item | Priority | Effort |
|---|---|---|
| Real-time EEG pattern recognition integration | Medium | 1-2 weeks |
| NeuroImaging AI (MONAI integration for lesion detection) | Low | 2-3 weeks |
| Natural language EHR note parsing for auto-scale scoring | Low | 1 week |
| Integration with HCLS AI Factory landing page | Low | 1 hour |

---

## 12. Relationship to HCLS AI Factory

The Neurology Intelligence Agent demonstrates the **time-critical clinical extension** of the HCLS AI Factory architecture. Neurological conditions — particularly acute stroke and status epilepticus — demand sub-minute decision support that combines imaging, clinical scales, genomics, and treatment guidelines simultaneously.

- **Same Milvus instance** — 13 new owned collections alongside existing `genomic_evidence` (read-only)
- **Same embedding model** — BGE-small-en-v1.5 (384-dim)
- **Same LLM** — Claude via Anthropic API
- **Same hardware** — NVIDIA DGX Spark ($4,699)
- **Same patterns** — Pydantic models, BaseIngestPipeline, knowledge graph, query expansion

The neurogenomics integration connects Stage 1 (genomic variants in epilepsy genes like SCN1A, or movement disorder genes like GBA/LRRK2) directly to this agent's clinical workflows for genotype-informed treatment selection.

---

## 13. Credits

- **Adam Jones**
- **Apache 2.0 License**

---

!!! warning "Clinical Decision Support Disclaimer"
    The Neurology Intelligence Agent is a clinical decision support research tool for neurological evaluation. It is not FDA-cleared and is not intended as a standalone diagnostic device. Time-critical decisions (stroke, status epilepticus) must follow institutional protocols. All recommendations should be reviewed by qualified neurologists. Apache 2.0 License.
