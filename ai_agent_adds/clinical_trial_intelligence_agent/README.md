# Clinical Trial Intelligence Agent


![Architecture Infographic](docs/images/infographic.jpg)

*Source: [clinical-trial-intelligence-agent](https://github.com/ajones1923/clinical-trial-intelligence-agent)*

RAG-powered clinical trial decision support agent for the HCLS AI Factory platform. Provides AI-driven protocol optimization, patient-trial matching, site selection, eligibility optimization, adaptive design evaluation, safety signal detection, regulatory document generation, competitive intelligence, diversity assessment, and decentralized trial planning.

## Architecture

```
                    +------------------+
                    |   Streamlit UI   |
                    |    (port 8128)   |
                    +--------+---------+
                             |
                    +--------v---------+
                    |   FastAPI Server  |
                    |    (port 8538)    |
                    +--------+---------+
                             |
              +--------------+--------------+
              |              |              |
     +--------v---+  +------v------+  +----v-------+
     | RAG Engine |  |  Workflow   |  |    LLM     |
     |            |  |   Engine    |  |  (Claude)  |
     +--------+---+  +------+------+  +----+-------+
              |              |              |
     +--------v--------------v--------------v------+
     |              Milvus Vector DB               |
     |           (14 collections, 384d)            |
     +---------------------------------------------+
```

## Knowledge Collections (14)

| Collection | Description |
|-----------|-------------|
| `trial_protocols` | Protocol designs, amendments, SOAs |
| `trial_eligibility` | Inclusion/exclusion criteria |
| `trial_endpoints` | Primary, secondary, exploratory endpoints |
| `trial_sites` | Site feasibility, enrollment history |
| `trial_investigators` | Investigator profiles, experience |
| `trial_results` | Published trial results, CSRs |
| `trial_regulatory` | FDA/EMA guidance, IND/NDA data |
| `trial_literature` | PubMed clinical trial publications |
| `trial_biomarkers` | Biomarker-driven trial designs |
| `trial_safety` | Safety data, DSMB reports, AE databases |
| `trial_rwe` | Real-world evidence, registries |
| `trial_adaptive` | Adaptive design references |
| `trial_guidelines` | ICH, FDA, EMA guidelines |
| `genomic_evidence` | Shared genomic evidence (cross-agent) |

## Workflows (10 + GENERAL)

1. **Protocol Design** -- Complexity scoring, SOA review, endpoint optimization
2. **Patient Matching** -- AI eligibility screening with genomic/biomarker integration
3. **Site Selection** -- Feasibility scoring, enrollment forecasting, diversity metrics
4. **Eligibility Optimization** -- Population impact modeling, competitor benchmarking
5. **Adaptive Design** -- Bayesian interim analysis, dose-response, futility assessment
6. **Safety Signal** -- Disproportionality analysis (PRR/ROR), causality assessment
7. **Regulatory Docs** -- IND, CSR, briefing doc generation (FDA/EMA/PMDA)
8. **Competitive Intel** -- Landscape analysis, enrollment race tracking
9. **Diversity Assessment** -- FDA diversity guidance compliance, gap analysis
10. **Decentralized Planning** -- DCT component feasibility, hybrid model design
11. **General** -- Free-form RAG Q&A across all collections

## Port Map

| Service | Port |
|---------|------|
| Streamlit UI | 8128 |
| FastAPI API | 8538 |
| Milvus (standalone) | 39530 |
| Milvus health | 39091 |

## Quickstart

### Local Development

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Configure environment
cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY

# 3. Start FastAPI server
uvicorn api.main:app --host 0.0.0.0 --port 8538 --reload

# 4. Start Streamlit UI (separate terminal)
streamlit run app/trial_ui.py --server.port 8128
```

### Docker Compose (Standalone)

```bash
cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY
docker compose up -d
docker compose logs -f trial-setup  # watch setup progress
```

### Verify

```bash
# Health check
curl http://localhost:8538/health

# List workflows
curl http://localhost:8538/workflows

# Query
curl -X POST http://localhost:8538/v1/trial/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What adaptive designs are used in Phase II oncology trials?"}'
```

## API Reference

See the auto-generated OpenAPI docs at `http://localhost:8538/docs` when the server is running.

## Tech Stack

- **LLM**: Claude (Anthropic) via `anthropic` SDK
- **Vector DB**: Milvus 2.4 with BGE-small-en-v1.5 embeddings (384d)
- **API**: FastAPI with Pydantic v2 schemas
- **UI**: Streamlit with NVIDIA dark theme
- **Container**: Multi-stage Docker, non-root `trialuser`
- **Monitoring**: Prometheus-compatible `/metrics` endpoint

## Author

Adam Jones -- HCLS AI Factory, March 2026

## License

This project is licensed under the [Apache License 2.0](LICENSE).
