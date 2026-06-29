# Neurology Intelligence Agent


![Architecture Infographic](docs/images/infographic.jpg)

*Source: [neurology-intelligence-agent](https://github.com/ajones1923/neurology-intelligence-agent)*

RAG-powered neurology clinical decision support system for the HCLS AI Factory.

## Features

- **8 Clinical Workflows**: Acute stroke triage, dementia evaluation, epilepsy classification, brain tumor grading, MS monitoring, Parkinson's assessment, headache classification, neuromuscular evaluation
- **10 Clinical Scale Calculators**: NIHSS, GCS, MoCA, MDS-UPDRS Part III, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn-Yahr
- **RAG-Powered Q&A**: Multi-collection semantic search across neurology knowledge base
- **Guideline Integration**: AAN, AHA/ASA, ILAE, ICHD-3, WHO CNS 2021, NCCN, McDonald 2017, MDS criteria
- **Multi-Format Reports**: Markdown, JSON, PDF, FHIR R4 DiagnosticReport
- **SSE Event Streaming**: Real-time workflow progress and cross-agent events

## Ports

| Service | Port |
|---------|------|
| FastAPI API | 8528 |
| Streamlit UI | 8529 |
| Milvus (standalone) | 59530 |

## Quick Start

```bash
# 1. Configure
cp .env.example .env
# Edit .env -- set ANTHROPIC_API_KEY

# 2. Start all services
docker compose up -d

# 3. Watch setup progress
docker compose logs -f neuro-setup

# 4. Access
# API:  http://localhost:8528/health
# UI:   http://localhost:8529
# Docs: http://localhost:8528/docs
```

## Local Development

```bash
# Install dependencies
pip install -r requirements.txt

# Start API
uvicorn api.main:app --host 0.0.0.0 --port 8528 --reload

# Start UI (separate terminal)
streamlit run app/neuro_ui.py --server.port 8529
```

## API Endpoints

### System
- `GET /health` -- Service health
- `GET /collections` -- Milvus collections
- `GET /workflows` -- Available workflows
- `GET /metrics` -- Prometheus metrics

### Clinical (prefix: `/v1/neuro`)
- `POST /query` -- RAG Q&A
- `POST /search` -- Multi-collection search
- `POST /scale/calculate` -- Clinical scale calculator
- `POST /stroke/triage` -- Acute stroke triage
- `POST /dementia/evaluate` -- Dementia evaluation
- `POST /epilepsy/classify` -- Epilepsy classification
- `POST /tumor/grade` -- Brain tumor grading
- `POST /ms/assess` -- MS assessment
- `POST /parkinsons/assess` -- Parkinson's assessment
- `POST /headache/classify` -- Headache classification
- `POST /neuromuscular/evaluate` -- Neuromuscular evaluation
- `POST /workflow/{type}` -- Generic workflow dispatch
- `GET /domains` -- Domain catalogue
- `GET /scales` -- Scale catalogue
- `GET /guidelines` -- Guideline reference
- `GET /knowledge-version` -- Version metadata

### Reports (`/v1/reports`)
- `POST /generate` -- Generate report
- `GET /formats` -- Supported formats

### Events (`/v1/events`)
- `GET /stream` -- SSE event stream
- `GET /health` -- SSE health

## Architecture

```
neurology_intelligence_agent/
  api/
    main.py              # FastAPI app (port 8528)
    routes/
      neuro_clinical.py  # Clinical endpoints + scale calculators
      reports.py         # Report generation
      events.py          # SSE streaming
  app/
    neuro_ui.py          # Streamlit UI (port 8529)
  config/
    settings.py          # Pydantic settings (NEURO_ env prefix)
  src/
    models.py            # Enums and data models
```

## Author

Adam Jones -- HCLS AI Factory, March 2026

## License

This project is licensed under the [Apache License 2.0](LICENSE).
