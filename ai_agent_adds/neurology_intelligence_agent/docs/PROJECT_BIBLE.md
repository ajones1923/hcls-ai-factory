# Neurology Intelligence Agent -- Project Bible

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones

---

## 1. Project Identity

**Name:** Neurology Intelligence Agent
**Codename:** NeuroAgent
**Repository:** `ai_agent_adds/neurology_intelligence_agent/`
**Parent Project:** HCLS AI Factory
**Platform Position:** Precision Intelligence Network -- one of three GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence Network, Therapeutic Discovery Engine)
**License:** Apache 2.0 (Open Source)

### Platform: 11 Intelligence Agents

| # | Agent | Domain |
|---|-------|--------|
| 1 | Precision Biomarker | Biomarker interpretation and trends |
| 2 | Precision Oncology | Cancer genomics and treatment |
| 3 | CAR-T Intelligence | CAR-T cell therapy |
| 4 | Medical Imaging | Multi-modal imaging analysis |
| 5 | Precision Autoimmune | Autoimmune disease analysis |
| 6 | Cardiology Intelligence | Cardiovascular decision support |
| 7 | Clinical Trial Intelligence | Trial matching and eligibility |
| 8 | **Neurology Intelligence** | **Neurological decision support** |
| 9 | Rare Disease Intelligence | Rare disease diagnosis |
| 10 | Pharmacogenomics (PGx) | Drug-gene interactions |
| 11 | Pediatric Oncology | Pediatric cancer intelligence |

### Mission Statement

Deliver an AI-powered clinical decision support system that unifies fragmented neurological evidence -- spanning cerebrovascular, neurodegenerative, epilepsy, movement disorders, multiple sclerosis, headache, neuromuscular, and neuro-oncology domains -- into a single RAG-driven intelligence platform. Clinicians receive guideline-grounded, evidence-cited recommendations in under five seconds.

### Why This Matters

- 3 billion people worldwide affected by neurological conditions (WHO 2024)
- Neurological data is fragmented across siloed subspecialties with inconsistent terminology
- Clinical decision-making requires synthesizing imaging, electrophysiology, genomics, and guideline evidence simultaneously
- Time-critical conditions (acute stroke, status epilepticus) demand sub-minute decision support

---

## 2. Architecture Summary

### Three-Layer Stack

```
Presentation     Streamlit UI (8529) + FastAPI REST (8528)
Intelligence     RAG Engine + Workflow Engine + Scale Calculators + Query Expansion
Data             Milvus 2.4 (14 collections, BGE-small 384-dim, IVF_FLAT/COSINE)
```

### Core Design Decisions

| Decision | Rationale |
|---|---|
| 14 domain-specific collections (not 1 monolithic) | Workflow-specific weight boosting yields higher-relevance results |
| BGE-small-en-v1.5 (384-dim) | Best size-accuracy tradeoff for medical text on DGX Spark |
| IVF_FLAT with nlist=128 | Sub-100ms search latency at 855K estimated records |
| Pydantic BaseSettings with `NEURO_` prefix | Type-safe config, env-driven, no secrets in code |
| ThreadPoolExecutor for parallel search | 14 collections searched concurrently; reduces latency from 14x to ~1x |
| Claude Sonnet for synthesis | Strongest medical reasoning among available LLMs |
| Hoehn-Yahr + MDS-UPDRS + MoCA + NIHSS + ... | 10 validated instruments prevent ad-hoc scoring |

---

## 3. Directory Structure

```
neurology_intelligence_agent/
  api/
    __init__.py
    main.py                       # FastAPI application with lifespan
    routes/
      __init__.py
      neuro_clinical.py           # 15 clinical endpoints
      reports.py                  # Report generation
      events.py                   # SSE streaming
  app/
    __init__.py
    neuro_ui.py                   # Streamlit chat interface
  config/
    __init__.py
    settings.py                   # Pydantic BaseSettings (50+ params)
  data/
    cache/                        # Conversation persistence
    reference/                    # Reference data files
  docs/
    NEUROLOGY_INTELLIGENCE_AGENT_RESEARCH_PAPER.md
    PRODUCTION_READINESS_REPORT.md
    PROJECT_BIBLE.md              # This document
    ARCHITECTURE_GUIDE.md
    WHITE_PAPER.md
    DEPLOYMENT_GUIDE.md
    DEMO_GUIDE.md
    LEARNING_GUIDE_FOUNDATIONS.md
    LEARNING_GUIDE_ADVANCED.md
    INDEX.md
  scripts/
    setup_collections.py          # Create Milvus schemas
    seed_knowledge.py             # Populate knowledge base
    run_ingest.py                 # Data ingestion pipeline
    generate_docx.py              # DOCX report generation
  src/
    __init__.py
    agent.py                      # Agent orchestrator
    clinical_scales.py            # 10 scale calculators
    clinical_workflows.py         # 8 clinical workflows
    collections.py                # 14 collection schemas
    cross_modal.py                # Cross-agent triggers
    export.py                     # Report formats
    knowledge.py                  # Domain knowledge base
    metrics.py                    # Prometheus metrics
    models.py                     # Enums and Pydantic models
    query_expansion.py            # 251+ aliases, 16 synonym maps
    rag_engine.py                 # Multi-collection RAG
    scheduler.py                  # Ingest scheduler
    utils/
      __init__.py
    ingest/
      __init__.py
      base.py                     # Base ingest pipeline
      pubmed_neuro_parser.py      # PubMed parser
      neuroimaging_parser.py      # Imaging protocol parser
      eeg_parser.py               # EEG pattern parser
  tests/
    __init__.py
    conftest.py                   # Shared fixtures
    test_agent.py                 # 5 tests
    test_api.py                   # 8 tests
    test_clinical_scales.py       # 35 tests
    test_clinical_workflows.py    # 11 tests
    test_collections.py           # 15 tests
    test_integration.py           # 16 tests
    test_knowledge.py             # 30 tests
    test_models.py                # 55 tests
    test_query_expansion.py       # 5 tests
    test_rag_engine.py            # 4 tests
    test_settings.py              # 18 tests
    test_workflow_execution.py    # 7 tests
  docker-compose.yml              # Standalone deployment
  Dockerfile                      # Container image
  README.md                       # Quick-start guide
  requirements.txt                # Python dependencies
```

---

## 4. Data Model

### Enums (18)

| Enum | Values | Purpose |
|---|---|---|
| `NeuroWorkflowType` | 9 values | Workflow routing |
| `NeuroDomain` | 10 values | Domain classification |
| `StrokeType` | 4 values | Stroke classification |
| `DementiaSubtype` | 9 values | Dementia differential |
| `ATNStage` | 8 values | Alzheimer's biomarker staging |
| `SeizureType` | 9 values | ILAE 2017 seizure types |
| `EpilepsySyndrome` | 10 values | Epilepsy syndrome identification |
| `MSPhenotype` | 4 values | MS clinical phenotypes |
| `DMTCategory` | 3 values | MS DMT efficacy tiers |
| `TumorGrade` | 4 values | WHO 2021 CNS grades |
| `TumorMolecularMarker` | 11 values | CNS tumor molecular markers |
| `ParkinsonsSubtype` | 3 values | PD motor subtypes |
| `HeadacheType` | 10 values | ICHD-3 headache types |
| `NMJPattern` | 8 values | EMG/NCS pattern classification |
| `SeverityLevel` | 5 values | Clinical severity |
| `EvidenceLevel` | 4 values | AAN evidence classification |
| `GuidelineClass` | 5 values | Recommendation classification |
| `ClinicalScaleType` | 10 values | Scale calculator dispatch |

### Pydantic Models (12)

`NeuroQuery`, `NeuroSearchResult`, `StrokeAssessment`, `DementiaAssessment`, `SeizureClassification`, `MSAssessment`, `TumorAssessment`, `HeadacheClassification`, `NeuromuscularAssessment`, `ScaleResult`, `WorkflowResult`, `NeuroResponse`

---

## 5. Key Technical Specifications

| Specification | Value |
|---|---|
| Embedding model | BAAI/bge-small-en-v1.5 |
| Embedding dimension | 384 |
| Vector index type | IVF_FLAT |
| Similarity metric | COSINE |
| IVF nlist | 128 |
| Score threshold | 0.4 |
| Max conversation context | 3 turns |
| Conversation TTL | 24 hours |
| Citation high threshold | 0.75 |
| Citation medium threshold | 0.60 |
| LLM model | claude-sonnet-4-6 |
| LLM max tokens | 2048 |
| LLM temperature | 0.7 |

---

## 6. Dependency Map

### Python Dependencies

| Package | Purpose |
|---|---|
| `fastapi` | REST API framework |
| `uvicorn` | ASGI server |
| `pydantic` + `pydantic-settings` | Data validation and configuration |
| `pymilvus` | Milvus vector database client |
| `sentence-transformers` | BGE-small embedding model |
| `anthropic` | Claude LLM client |
| `streamlit` | Chat UI |
| `loguru` | Structured logging |
| `python-docx` | DOCX report generation |

### Infrastructure Dependencies

| Service | Version | Purpose |
|---|---|---|
| Milvus | 2.4 | Vector database |
| etcd | 3.5.5 | Milvus metadata |
| MinIO | 2023.03 | Milvus object storage |
| Docker | 24+ | Container runtime |
| Python | 3.10+ | Runtime |

---

## 7. Naming Conventions

| Element | Convention | Example |
|---|---|---|
| Collection names | `neuro_` prefix, snake_case | `neuro_cerebrovascular` |
| Environment variables | `NEURO_` prefix, UPPER_SNAKE | `NEURO_MILVUS_HOST` |
| Enum values | snake_case strings | `acute_stroke` |
| API paths | lowercase, kebab-compatible | `/v1/neuro/stroke/triage` |
| Test files | `test_` prefix | `test_clinical_scales.py` |
| Scale classes | PascalCase + Calculator suffix | `NIHSSCalculator` |
| Workflow classes | PascalCase + Workflow suffix | `AcuteStrokeWorkflow` |

---

## 8. Guideline Bodies Referenced

| Organization | Abbreviation | Domains |
|---|---|---|
| American Academy of Neurology | AAN | All neurology domains |
| European Academy of Neurology | EAN | Cross-domain European guidelines |
| American Heart/Stroke Association | AHA/ASA | Cerebrovascular disease |
| International League Against Epilepsy | ILAE | Epilepsy classification and treatment |
| International Headache Society | IHS | ICHD-3 headache classification |
| Movement Disorder Society | MDS | PD, dystonia, tremor |
| National Comprehensive Cancer Network | NCCN | CNS tumors |
| World Health Organization | WHO | CNS tumor classification (2021) |
| NIA-AA | NIA-AA | ATN framework for Alzheimer's |
| American Clinical Neurophysiology Society | ACNS | EEG terminology |

---

## 9. Quality Gates

| Gate | Criterion | Status |
|---|---|---|
| Unit tests pass | 209/209 tests green | Pass |
| Scale calculator accuracy | All 10 calculators match clinical definitions | Pass |
| Collection weight sums | Each workflow sums to ~1.0 (tolerance 0.02) | Pass |
| Model serialization | All 12 Pydantic models serialize/deserialize cleanly | Pass |
| API contract | All endpoints return documented response schemas | Pass |
| Configuration validation | Settings.validate() returns no critical issues | Pass |
| Docker build | `docker compose build` succeeds | Pass |
| Knowledge integrity | Drug count >= 42, gene count >= 38 | Pass |

---

## 10. Cross-Agent Integration and Pediatric Oncology Focus

### Cross-Agent Endpoints

The Neurology Intelligence Agent calls 5 sibling agents and exposes a `/v1/neuro/integrated-assessment` endpoint:

- **Imaging Agent (8524)**: MRI lesion characterization, DWI/FLAIR mismatch, tumor segmentation
- **Cardiology Agent (8126)**: Cardioembolic stroke evaluation (AF detection, PFO, echocardiography)
- **Biomarker Agent (8529)**: CSF biomarker trending (amyloid, tau, NfL trajectories)
- **Clinical Trial Agent (8537)**: Trial matching for anti-amyloid, gene therapy, novel DMTs
- **Rare Disease Agent (8134)**: Rare neurogenetic disorder evaluation (inherited ataxias, leukodystrophies)

### Pediatric Neuro-Oncology Focus

The neurology agent provides specialized support for neurotoxicity in pediatric oncology:

- **MTX leukoencephalopathy**: 3-10% of pediatric ALL patients; white matter changes on MRI correlated with intrathecal/high-dose MTX exposure
- **Vincristine neuropathy**: 30-40% incidence; peripheral neuropathy tracking with dose modification recommendations
- **PRES**: Posterior reversible encephalopathy syndrome detection in pediatric transplant patients
- **L-asparaginase thrombosis**: CNS thrombotic event monitoring (sagittal sinus, cortical vein thrombosis)
- **Cranial radiation effects**: Long-term neurocognitive monitoring, leukoencephalopathy, secondary tumor risk
- **Pediatric CNS tumors**: Medulloblastoma molecular subgrouping (WNT/SHH/Group 3/Group 4) and diffuse midline glioma (H3K27M-altered) support

---

## 11. Roadmap

### v1.0 (Current -- March 2026)

- 14 collections, 8 workflows, 10 scales
- Standalone and integrated Docker deployment
- 209 tests across 12 modules
- Full documentation set

### v1.1 (Planned -- Q2 2026)

- EHR integration via FHIR R4 DiagnosticReport
- Real-time EEG pattern classification
- DWI/FLAIR MRI feature extraction
- Structured medication interaction checking

### v1.2 (Planned -- Q3 2026)

- Multi-institutional knowledge federation
- Longitudinal patient tracking
- Automated guideline update ingestion
- Clinical note summarization

---

## 12. Contact and Ownership

| Role | Contact |
|---|---|
| Project Lead / Developer | Adam Jones |
| Platform | HCLS AI Factory |
| Repository | `ai_agent_adds/neurology_intelligence_agent/` |
| Documentation | `docs/` directory |
| Issues | GitHub Issues on hcls-ai-factory |

---

*Neurology Intelligence Agent -- Project Bible v1.3.0*
*HCLS AI Factory / GTC Europe 2026*
