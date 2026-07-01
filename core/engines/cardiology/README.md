# Cardiology Intelligence Agent


![Architecture Infographic](docs/images/infographic.jpg)

*Source: [cardiology-intelligence-agent](https://github.com/ajones1923/cardiology-intelligence-agent)*

RAG-powered cardiovascular clinical decision support system built on Milvus, Claude, and BGE-small-en-v1.5. Part of the HCLS AI Factory precision medicine platform.

Synthesizes cardiac imaging, electrophysiology, hemodynamics, heart failure management, valvular disease, preventive cardiology, interventional data, and cardio-oncology surveillance into guideline-aligned clinical recommendations using ACC/AHA/ESC evidence.

**Author:** Adam Jones
**Date:** March 2026

## Architecture

```
                    +------------------+
                    |  Streamlit UI    |
                    |  :8536           |
                    +--------+---------+
                             |
                    +--------+---------+
                    |  FastAPI Server  |
                    |  :8126           |
                    +--------+---------+
                             |
              +--------------+--------------+
              |                             |
    +---------+----------+     +------------+-----------+
    |  RAG Engine        |     |  Clinical Engines      |
    |  - Query Expansion |     |  - Risk Calculators    |
    |  - Multi-Collection|     |  - GDMT Optimizer      |
    |  - Citation Scoring|     |  - Cross-Modal Triggers|
    +---------+----------+     +------------+-----------+
              |                             |
    +---------+-----------------------------+-----------+
    |              Milvus Vector DB (:19530)            |
    |  12 collections + genomic_evidence (384-dim)     |
    +--------------------------------------------------+
    |  etcd (:2379)  |  MinIO (:9000)                  |
    +--------------------------------------------------+
```

## Collections (13)

| Collection | Description |
|---|---|
| cardio_literature | Published cardiovascular research, reviews, meta-analyses |
| cardio_trials | Cardiovascular clinical trials and landmark trial results |
| cardio_imaging | Cardiac imaging protocols, findings, measurements (echo, CT, MRI, nuclear) |
| cardio_electrophysiology | ECG interpretation, arrhythmia classification, EP data |
| cardio_heart_failure | HF classification, GDMT protocols, management algorithms |
| cardio_valvular | Valvular heart disease assessment and intervention criteria |
| cardio_prevention | Preventive cardiology: risk stratification, lipid management |
| cardio_interventional | Interventional procedures, techniques, outcomes |
| cardio_oncology | Cardio-oncology surveillance, cardiotoxicity detection |
| cardio_devices | FDA-cleared cardiovascular AI devices, implantable devices |
| cardio_guidelines | ACC/AHA/ESC/HRS clinical practice guidelines |
| cardio_hemodynamics | Catheterization data, pressure tracings, derived calculations |
| genomic_evidence | Shared genomic evidence (read-only, 3.5M variants) |

## Port Map

| Service | Port |
|---|---|
| Streamlit UI | 8536 |
| FastAPI API | 8126 |
| Milvus gRPC | 19530 |
| Milvus Health | 9091 |
| MinIO | 9000 |
| etcd | 2379 |

## Quickstart

```bash
# 1. Configure environment
cp .env.example .env
# Edit .env and set ANTHROPIC_API_KEY

# 2. Start all services (Milvus + Cardio agent)
docker compose up -d

# 3. Watch setup/seed progress
docker compose logs -f cardio-setup

# 4. Open the UI
open http://localhost:8536
```

### Manual Setup (without Docker)

```bash
# Install dependencies
pip install -r requirements.txt

# Create Milvus collections
python scripts/setup_collections.py --drop-existing --seed

# Seed knowledge base
python scripts/seed_knowledge.py

# Run live data ingest (PubMed, ClinicalTrials.gov, ACC/AHA, etc.)
python scripts/run_ingest.py

# Start Streamlit UI
streamlit run app/cardio_ui.py --server.port=8536

# Start FastAPI server (separate terminal)
uvicorn api.main:app --host 0.0.0.0 --port 8126 --workers 2
```

## Clinical Workflows (11)

1. **Coronary Artery Disease Assessment** -- Calcium scoring, CAD-RADS, plaque characterization
2. **Heart Failure Classification & GDMT** -- NYHA/ACC staging, 4-pillar GDMT optimization
3. **Valvular Heart Disease Quantification** -- Severity grading, intervention criteria
4. **Arrhythmia Detection & Management** -- ECG interpretation, CHA2DS2-VASc, anticoagulation
5. **Cardiac MRI Tissue Characterization** -- LGE patterns, T1/T2 mapping, parametric analysis
6. **Stress Test Interpretation** -- Duke Treadmill Score, SPECT/PET MPI, stress echo
7. **Preventive Risk Stratification** -- ASCVD risk, statin eligibility, CAC reclassification
8. **Cardio-Oncology Surveillance** -- Cardiotoxicity monitoring, CTRCD detection, GLS tracking
9. **Acute Decompensated HF** -- Acute HF triage, diuretic dosing, hemodynamic profiling (warm/cold, wet/dry)
10. **Post-MI** -- Post-MI risk stratification, secondary prevention, cardiac rehabilitation planning
11. **Myocarditis/Pericarditis** -- Inflammatory cardiac disease assessment, Lake Louise criteria, treatment protocols

## Risk Calculators (6)

- **ASCVD (PCE)** -- 10-year atherosclerotic cardiovascular disease risk
- **HEART Score** -- Chest pain risk stratification in ED
- **CHA2DS2-VASc** -- Atrial fibrillation stroke risk
- **HAS-BLED** -- Anticoagulation bleeding risk
- **MAGGIC** -- Heart failure mortality risk
- **EuroSCORE II** -- Cardiac surgical mortality risk

## Data Sources

- **ACC/AHA** -- American College of Cardiology / American Heart Association guidelines
- **ESC** -- European Society of Cardiology guidelines
- **PubMed** -- Cardiovascular literature via NCBI E-utilities
- **ClinicalTrials.gov** -- Cardiovascular clinical trials
- **ASE/SCCT/SCMR** -- Imaging society guideline documents
- **FDA** -- FDA-cleared cardiovascular AI/ML devices

## Tech Stack

- **Vector DB:** Milvus 2.4 (IVF_FLAT / COSINE)
- **Embeddings:** BGE-small-en-v1.5 (384-dim)
- **LLM:** Claude Sonnet 4.6 (Anthropic)
- **UI:** Streamlit
- **API:** FastAPI + Uvicorn
- **Monitoring:** Prometheus metrics
- **Compute:** NVIDIA DGX Spark

## License

This project is licensed under the [Apache License 2.0](LICENSE).

## Author

Adam Jones -- HCLS AI Factory, March 2026
