# Biomarker Intelligence Agent


![Architecture Infographic](docs/images/infographic.jpg)

*Source: [precision-biomarker-agent](https://github.com/ajones1923/precision-biomarker-agent)*

**HCLS AI Factory** | Genotype-aware biomarker interpretation, biological age estimation, disease trajectory detection, and pharmacogenomic profiling.

## Quick Start

```bash
# 1. Configure API key
cp .env.example .env
# Edit .env and add ANTHROPIC_API_KEY=sk-ant-...

# 2. Launch all services (Milvus + Streamlit + API + seed)
docker compose up -d

# 3. Watch seed progress
docker compose logs -f biomarker-setup

# 4. Open the UI
open http://localhost:8528
```

### Local Development (no Docker)

```bash
pip install -r requirements.txt

# Start Milvus (must be running on localhost:19530)
python scripts/setup_collections.py --drop-existing
python scripts/seed_all.py
python scripts/validate_e2e.py

# Launch Streamlit UI
streamlit run app/biomarker_ui.py --server.port 8528
```

## Architecture Overview

```
Patient Biomarkers + Genotypes
        |
        v
+-------------------+     +---------------------+
|  Biological Age   |     |  Disease Trajectory  |
|  (PhenoAge/Grim)  |     |  (6 disease models)  |
+-------------------+     +---------------------+
        |                          |
        v                          v
+-------------------+     +---------------------+
| Pharmacogenomics  |     | Genotype Adjustment  |
| (7 PGx genes)     |     | (reference ranges)   |
+-------------------+     +---------------------+
        |                          |
        +----------+  +----------+
                   |  |
                   v  v
           +------------------+
           |   RAG Engine     |
           | (10 collections) |
           +------------------+
                   |
                   v
           +------------------+
           | Report Generator |
           | (PDF / FHIR R4)  |
           +------------------+
```

## Milvus Collections (10 Domain + 1 Shared)

| # | Collection                       | Description                                      | Records |
|---|----------------------------------|--------------------------------------------------|---------|
| 1 | `biomarker_reference`            | Reference biomarker definitions and ranges        | ~60     |
| 2 | `biomarker_genetic_variants`     | Genetic variants affecting biomarker levels       | ~30     |
| 3 | `biomarker_pgx_rules`            | CPIC pharmacogenomic dosing rules                 | ~50     |
| 4 | `biomarker_disease_trajectories` | Disease progression stage definitions             | ~30     |
| 5 | `biomarker_clinical_evidence`    | Published clinical evidence                       | ~40     |
| 6 | `biomarker_nutrition`            | Genotype-aware nutrition guidelines               | ~25     |
| 7 | `biomarker_drug_interactions`    | Gene-drug interactions                            | ~35     |
| 8 | `biomarker_aging_markers`        | Epigenetic aging clock markers (PhenoAge/GrimAge) | ~15     |
| 9 | `biomarker_genotype_adjustments` | Genotype-based reference range adjustments        | ~20     |
|10 | `biomarker_monitoring`           | Condition-specific monitoring protocols           | ~15     |
|+1 | `genomic_evidence`               | Shared genomic variants (read-only)               | varies  |

All collections use BGE-small-en-v1.5 (384-dim) embeddings with IVF_FLAT/COSINE indexing.

## Analysis Modules

1. **Biological Age Engine** -- PhenoAge (Levine 2018) and GrimAge surrogate estimation from 9 routine blood biomarkers.
2. **Disease Trajectory Analyzer** -- Pre-symptomatic detection across 6 disease categories with genotype-stratified thresholds.
3. **Pharmacogenomic Mapper** -- Star allele interpretation for 7 pharmacogenes (CYP2D6, CYP2C19, CYP2C9, CYP3A5, SLCO1B1, VKORC1, MTHFR) plus HLA-B*57:01.
4. **Genotype Adjustment Engine** -- Adjusts standard reference ranges based on patient genotype (PNPLA3, TCF7L2, APOE, etc.).
5. **RAG Evidence Engine** -- Cross-collection semantic search with Claude-powered synthesis.
6. **Nutrition Advisor** -- Genotype-aware supplementation recommendations.
7. **Report Generator** -- 12-section clinical report with PDF and FHIR R4 export.

## API Endpoints

| Method | Path                  | Description                           |
|--------|-----------------------|---------------------------------------|
| GET    | `/health`             | Service health check                  |
| POST   | `/analyze`            | Full patient analysis                 |
| POST   | `/biological-age`     | PhenoAge/GrimAge calculation          |
| POST   | `/disease-risk`       | Disease trajectory analysis           |
| POST   | `/pgx`                | Pharmacogenomic profile               |
| POST   | `/search`             | RAG evidence search                   |
| POST   | `/report`             | Generate clinical report              |
| GET    | `/collections/stats`  | Collection record counts              |

## Port Assignments

| Port | Service                   |
|------|---------------------------|
| 8528 | Streamlit UI              |
| 8529 | FastAPI REST API           |
| 19530| Milvus gRPC               |
| 9091 | Milvus health/metrics     |

## Tests

```bash
pytest tests/ -v
```

All tests are self-contained and use mocks (no Milvus or LLM required).

## License

Apache 2.0
