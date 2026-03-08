# Precision Biomarker Intelligence Agent

Genomics-informed biomarker interpretation with 9 clinical analysis engines. Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory).

## Overview

The Precision Biomarker Intelligence Agent transforms raw lab results into clinically actionable intelligence by combining multi-collection RAG search with deterministic clinical engines. It integrates patient genotype data (ApoE, MTHFR, CYP variants), pharmacogenomic star alleles, and age/sex-stratified reference ranges to produce personalized interpretations that a standard lab report cannot provide.

| Collection | Records | Content |
|---|---|---|
| **Biomarker Reference** | 208 | Comprehensive biomarker definitions, units, and clinical significance |
| **Genetic Variants** | 42 | Clinically actionable SNPs with risk alleles and effect sizes |
| **PGx Rules** | 29 | CPIC Level 1A pharmacogenomic guidelines |
| **Disease Trajectories** | 39 | Multi-biomarker risk patterns for 9 disease domains |
| **Clinical Evidence** | 80 | Published clinical study evidence for biomarker interpretation |
| **Nutrition** | 50 | Nutrient-biomarker interactions and dietary recommendations |
| **Drug Interactions** | 51 | Medication effects on biomarker levels |
| **Aging Markers** | 20 | PhenoAge and GrimAge epigenetic clock biomarkers |
| **Genotype Adjustments** | 30 | Genotype-specific reference range modifications |
| **Monitoring** | 30 | Follow-up testing schedules and monitoring protocols |
| **Critical Values** | 21 | Threshold rules for critical/urgent/warning alerts |
| **Discordance Rules** | 12 | Cross-biomarker discordance detection patterns |
| **AJ Carrier Screening** | 10 | Ashkenazi Jewish genetic carrier screening panel |
| **Genomic Evidence** | 30 | *(read-only)* Shared from Stage 2 RAG pipeline |
| **Total** | **652 vectors** | **14 collections (13 owned + 1 read-only)** |

### 9 Clinical Analysis Engines

| Engine | Function |
|---|---|
| **BiologicalAgeCalculator** | PhenoAge (Levine 2018) and GrimAge (Lu 2019) composite scoring with confidence intervals |
| **DiseaseTrajectoryAnalyzer** | 9-domain risk stratification (cardiovascular, diabetes, liver, thyroid, iron, nutritional, kidney, bone, cognitive) |
| **PharmacogenomicMapper** | CPIC Level 1A star allele and genotype-to-phenotype mapping for 13 genes |
| **GenotypeAdjuster** | Genotype-aware reference range modifications + age-stratified adjustments (5 brackets) |
| **CriticalValueEngine** | Real-time critical/urgent/warning threshold detection with escalation routing |
| **DiscordanceDetector** | Cross-biomarker pattern analysis for clinically discordant results |
| **LabRangeInterpreter** | Three-way comparison: Quest vs LabCorp vs Function Health optimal ranges |
| **AJ Carrier Screening** | Ashkenazi Jewish population-specific genetic carrier screening (via RAG) |
| **Age-Stratified Adjustments** | Sex-stratified reference ranges across 5 age brackets (0-17, 18-39, 40-59, 60-79, 80+) |

### Example Queries

```
"Interpret my LDL of 138 given ApoE E3/E4 genotype"
"What does MTHFR C677T heterozygous mean for my homocysteine levels?"
"CYP2D6 *1/*4 — which medications should I avoid?"
"How does my HbA1c of 5.6 with TCF7L2 CT genotype affect diabetes risk?"
"My ferritin is 28 — is this normal for a 38-year-old female?"
```

### Demo Guide

For a complete walkthrough of all 8 UI tabs with both sample patients, see the **[Demo Guide](demo-guide.md)**.

## Architecture

```
Patient Data Input
    |
    v
[Critical Value Engine] ──── CRITICAL? ──── YES ──> Immediate Alert
    |                                                (escalation routing)
    NO
    |
    v
[9 Clinical Analysis Engines — Parallel Execution]
    |               |              |              |
    v               v              v              v
BiologicalAge   Disease       Pharmaco-      Genotype
Calculator      Trajectory    genomic        Adjuster
(PhenoAge +     Analyzer      Mapper         (age-stratified
 GrimAge)       (9 domains)   (13 genes)      ranges)
    |               |              |              |
    +-------+-------+--------------+--------------+
            |
            v
    [Multi-Collection RAG Engine]
    Parallel search across 14 Milvus collections
    (ThreadPoolExecutor, configurable weights)
            |
            v
    [Claude Sonnet 4.6] -> Grounded response with citations
            |
            v
    [Export: FHIR R4 | PDF | Markdown | CSV]
```

Built on the HCLS AI Factory platform:

- **Vector DB:** Milvus 2.4 with IVF_FLAT/COSINE indexes (nlist=1024, nprobe=16)
- **Embeddings:** BGE-small-en-v1.5 (384-dim)
- **LLM:** Claude Sonnet 4.6 (Anthropic API)
- **UI:** Streamlit 8 tabs (port 8528) | **API:** FastAPI (port 8529)
- **Hardware target:** NVIDIA DGX Spark ($3,999)

### UI Tabs

| # | Tab | Content |
|---|---|---|
| 1 | Biomarker Analysis | Critical alerts, discordance, lab range comparison, full analysis |
| 2 | Biological Age | PhenoAge + GrimAge with aging drivers and confidence intervals |
| 3 | Disease Risk | 9-domain risk cards with genotype integration |
| 4 | PGx Profile | Star allele and genotype-to-drug mapping (13 genes) |
| 5 | Evidence Explorer | RAG search across 14 collections with collection filtering |
| 6 | Reports | FHIR R4, PDF, Markdown, CSV export with validation |
| 7 | Patient 360 | Unified cross-agent dashboard (genomics + biomarkers + drugs) |
| 8 | Longitudinal | Multi-visit biomarker trending and trajectory analysis |

## Setup

### Prerequisites

- Python 3.10+
- Milvus 2.4 running on `localhost:19530`
- `BIOMARKER_ANTHROPIC_API_KEY` environment variable (or in `.env`)

### Install

```bash
cd ai_agent_adds/precision_biomarker_agent
pip install -r requirements.txt
```

### 1. Create Collections and Seed Reference Data

```bash
python3 scripts/seed_all.py
```

Creates 14 Milvus collections with IVF_FLAT indexes and seeds 652 vectors from 14 JSON reference files.

### 2. Validate Data Layer

```bash
python3 scripts/validate_e2e.py
```

Tests collection existence, vector counts, embedding search, and cross-collection queries.

### 3. Run Demo Validation

```bash
python3 scripts/demo_validation.py
```

Runs 65 checks across all 8 UI tabs: critical values, discordance detection, lab ranges, biological age (PhenoAge/GrimAge), disease trajectories (9 domains), pharmacogenomics (CYP2D6/CYP2C19), FHIR export with validation, genotype adjustments, and sample patient data integrity.

### 4. Run Unit Tests

```bash
python3 -m pytest tests/ -v
```

709 tests covering all modules: critical values (29), discordance (25), lab ranges (36), collections (125), UI (44), RAG engine, edge cases, and more.

### 5. Launch UI

```bash
streamlit run app/biomarker_ui.py --server.port 8528
```

### 6. Launch API (separate terminal)

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8529
```

## Project Structure

```
precision_biomarker_agent/
├── src/
│   ├── models.py                  # Pydantic data models (20+ models)
│   ├── collections.py             # 14 Milvus collection schemas + manager
│   ├── rag_engine.py              # Multi-collection RAG engine + Claude
│   ├── agent.py                   # PrecisionBiomarkerAgent orchestrator
│   ├── biological_age.py          # PhenoAge + GrimAge calculators
│   ├── disease_trajectory.py      # 9-domain risk stratification
│   ├── pharmacogenomics.py        # CPIC PGx mapper (13 genes)
│   ├── genotype_adjustment.py     # Genotype + age-stratified adjustments
│   ├── critical_values.py         # Critical/urgent/warning threshold engine
│   ├── discordance_detector.py    # Cross-biomarker discordance detection
│   ├── lab_range_interpreter.py   # Quest vs LabCorp vs optimal ranges
│   └── export.py                  # FHIR R4, PDF, Markdown, CSV export
├── app/
│   └── biomarker_ui.py            # Streamlit UI (8 tabs, ~1,700 lines)
├── api/
│   └── main.py                    # FastAPI REST server
├── config/
│   └── settings.py                # Pydantic BaseSettings (14 collection weights)
├── data/
│   └── reference/
│       ├── biomarker_reference.json           # 208 biomarker definitions
│       ├── biomarker_genetic_variants.json     # 42 clinically actionable SNPs
│       ├── biomarker_pgx_rules.json           # 29 CPIC pharmacogenomic rules
│       ├── biomarker_disease_trajectories.json # 39 disease trajectory patterns
│       ├── biomarker_clinical_evidence.json    # 80 clinical evidence records
│       ├── biomarker_nutrition.json            # 50 nutrient-biomarker interactions
│       ├── biomarker_drug_interactions.json    # 51 drug-biomarker effects
│       ├── biomarker_aging_markers.json        # 20 epigenetic clock markers
│       ├── biomarker_genotype_adjustments.json # 30 genotype-specific adjustments
│       ├── biomarker_monitoring.json           # 30 follow-up protocols
│       ├── biomarker_critical_values.json      # 21 critical value rules
│       ├── biomarker_discordance_rules.json    # 12 discordance patterns
│       ├── biomarker_aj_carrier_screening.json # 10 AJ carrier screening records
│       ├── biomarker_genomic_evidence.json     # 30 genomic evidence records
│       ├── biomarker_lab_ranges.json           # Multi-lab reference ranges
│       └── biomarker_sample_patients.json      # 2 sample patients (Male 45, Female 38)
├── scripts/
│   ├── seed_all.py                # Seed all 14 Milvus collections
│   ├── validate_e2e.py            # End-to-end data layer validation
│   └── demo_validation.py         # 65-check demo validation script
├── tests/
│   ├── conftest.py                # Shared fixtures (mock collections, patients)
│   ├── test_critical_values.py    # 29 tests
│   ├── test_discordance_detector.py # 25 tests
│   ├── test_lab_range_interpreter.py # 36 tests
│   ├── test_collections.py        # 125 tests
│   ├── test_ui.py                 # 44 tests
│   ├── test_rag_engine.py         # RAG engine tests
│   └── test_edge_cases.py         # Edge case coverage
├── docker-compose.yml
├── Dockerfile
├── requirements.txt
└── LICENSE                        # Apache 2.0
```

**57 Python files | ~29,000 lines | Apache 2.0**

## Performance

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified memory):

| Metric | Value |
|---|---|
| Seed all 14 collections (652 vectors) | ~2 min |
| Vector search (14 collections, top-5 each) | 10-18 ms (cached) |
| PhenoAge + GrimAge calculation | <100 ms |
| Disease trajectory analysis (9 domains) | <200 ms |
| Pharmacogenomic mapping (all genes) | <50 ms |
| Critical value check (21 rules) | <10 ms |
| Full RAG query (search + Claude) | ~20-30 sec |
| FHIR R4 export + validation | <500 ms |
| Demo validation (65 checks) | ~15 sec |
| Unit tests (709 tests) | ~45 sec |

## Status

- **Phase 1 (Scaffold)** -- Complete. Architecture, data models, 14 collection schemas, 9 clinical engines, RAG engine, agent orchestrator, and Streamlit UI.
- **Phase 2 (Data)** -- Complete. 652 vectors seeded across 14 Milvus collections from curated reference JSON files. 2 sample patients with full biomarkers, genotypes, star alleles, and clinical context.
- **Phase 3 (Integration)** -- Complete. Full RAG pipeline with Claude generating grounded, genotype-aware interpretations. FHIR R4 export with structural validation.
- **Phase 4 (Testing)** -- Complete. 709 unit tests passing. 65-check demo validation covering all 8 UI tabs. End-to-end data layer validation.
- **Phase 5 (Demo Ready)** -- Complete. Production-quality Streamlit UI with 8 tabs, sample patient auto-load, critical alert banners, and multi-format export.

## Credits

- **Adam Jones**
- **Apache 2.0 License**
