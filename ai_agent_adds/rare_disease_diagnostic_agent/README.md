# Rare Disease Diagnostic Agent


![Architecture Infographic](docs/images/infographic.jpg)

*Source: [rare-disease-diagnostic-agent](https://github.com/ajones1923/rare-disease-diagnostic-agent)*

RAG-powered rare disease diagnostic intelligence agent for the HCLS AI Factory.
Provides differential diagnosis, ACMG/AMP variant interpretation, HPO-based
phenotype matching, therapeutic option search, and clinical trial eligibility
across 14 specialized vector collections.

## Architecture

```
Patient Phenotype/Genotype
        |
        v
  [FastAPI REST API :8134]
        |
        +-- Differential Diagnosis (HPO + OMIM + Orphanet)
        +-- Variant Interpretation (ACMG/AMP + ClinVar)
        +-- Phenotype Matching (IC-weighted similarity)
        +-- Therapeutic Search (orphan drugs, gene therapy, ERT)
        +-- Trial Eligibility (ClinicalTrials.gov)
        +-- Case Similarity (literature + case reports)
        |
        v
  [Milvus Vector DB :19530]     [Claude LLM]
  14 collections                 Evidence synthesis
        |
        v
  [Streamlit UI :8544]
  5-tab diagnostic interface
```

## Quick Start

```bash
# 1. Configure
cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY

# 2. Start (standalone with Milvus)
docker compose up -d

# 3. Access
# API:  http://localhost:8134/health
# UI:   http://localhost:8544
# Docs: http://localhost:8134/docs
```

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with component status |
| GET | `/collections` | Loaded collection names |
| GET | `/workflows` | Available diagnostic workflows (11) |
| GET | `/metrics` | Prometheus-compatible metrics |
| POST | `/v1/diagnostic/query` | RAG Q&A query |
| POST | `/v1/diagnostic/search` | Multi-collection search |
| POST | `/v1/diagnostic/diagnose` | Differential diagnosis |
| POST | `/v1/diagnostic/variants/interpret` | ACMG variant classification |
| POST | `/v1/diagnostic/phenotype/match` | HPO-to-disease matching |
| POST | `/v1/diagnostic/therapy/search` | Therapeutic option search |
| POST | `/v1/diagnostic/trial/match` | Clinical trial eligibility |
| POST | `/v1/diagnostic/workflow/{type}` | Generic workflow dispatch |
| GET | `/v1/diagnostic/disease-categories` | Disease category catalog |
| GET | `/v1/diagnostic/gene-therapies` | Approved gene therapies |
| GET | `/v1/diagnostic/acmg-criteria` | ACMG criteria reference |
| GET | `/v1/diagnostic/hpo-categories` | HPO top-level terms |
| GET | `/v1/diagnostic/knowledge-version` | Knowledge base metadata |
| POST | `/v1/reports/generate` | Report generation |
| GET | `/v1/reports/formats` | Supported export formats |
| GET | `/v1/events/stream` | SSE event stream |
| GET | `/v1/events/health` | SSE subsystem health |

## Streamlit UI Tabs

1. **Patient Intake** -- HPO term entry, clinical notes, VCF upload
2. **Diagnostic Dashboard** -- Differential diagnosis with confidence scores
3. **Variant Review** -- ACMG/AMP classification with criteria breakdown
4. **Therapeutic Options** -- Approved therapies, gene therapy, trials
5. **Report Generator** -- Export as Markdown, JSON, PDF, FHIR, Phenopacket

## Knowledge Collections

14 specialized Milvus collections prefixed with `rd_`:
phenotypes, diseases, genes, variants, literature, trials, therapies,
case_reports, guidelines, pathways, registries, natural_history,
newborn_screening, plus shared `genomic_evidence`.

## Ports

| Service | Port |
|---------|------|
| FastAPI API | 8134 |
| Streamlit UI | 8544 |
| Milvus (standalone) | 49530 |
| Milvus metrics | 49091 |

## Development

```bash
# API server (hot reload)
uvicorn api.main:app --host 0.0.0.0 --port 8134 --reload

# Streamlit UI
streamlit run app/diagnostic_ui.py --server.port 8544

# Verify import
python3 -c "from api.main import app; print(app.title)"
```

## Report Formats

- **Markdown** -- Human-readable diagnostic report
- **JSON** -- Structured data export
- **PDF** -- Printable clinical report (placeholder)
- **FHIR R4** -- HL7 FHIR DiagnosticReport resource
- **GA4GH Phenopacket** -- Phenopacket v2 for rare disease data exchange

## Author

Adam Jones -- HCLS AI Factory, March 2026

## License

This project is licensed under the [Apache License 2.0](LICENSE).
