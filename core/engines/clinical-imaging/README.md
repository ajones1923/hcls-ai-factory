# Clinical Imaging Engine

Engine 4 of the HCLS AI Factory. Live DICOM analysis with cross-modality reasoning over images
and genomics. Default port **8525**.

## What it does
- **Segmentation & classification** — VISTA-3D, MAISI, and VILA-M3 workflows over CT / MRI / CXR.
- **Cross-modal reasoning** — joins imaging findings with the genomic and biomarker layers.
- **FHIR R4 export** — structured, interoperable results.

## Layout
- **`agent/`** — the deployable imaging application (FastAPI + UI, DICOM ingest, model serving,
  its own `docker-compose.yml`). This is the runnable service. Its heavy model weights and DICOM
  data live under `agent/data/` and are gitignored (hydrate locally).
- **`docs/`** — design, deployment, and demo guides for the engine.
- **`generate_imaging_report.py`** — standalone clinical imaging report generator.

## Run
See `agent/README.md` for the full service bring-up. In brief:
```bash
cd agent
docker compose up -d      # brings up the imaging service on :8525
```

## Roadmap (this engine)
Broader modality coverage, tighter genomics ↔ imaging joins, and registration of each imaging
workflow as a governed capability in the platform registry.
