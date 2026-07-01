# Precision Oncology Engine

Engine 5 of the HCLS AI Factory. Molecular Tumor Board (MTB) packet generation, therapy ranking,
and trial matching, with cross-modal joins to the imaging, biomarker, and trial capabilities.
Default port **8526**.

## What it does
- **MTB packet generation** — assemble a molecular tumor board packet from variants + evidence.
- **Therapy ranking** — rank therapies against the patient's molecular profile.
- **Trial matching** — match patients to clinical trials.
- **Cross-modal joins** — pull in imaging, biomarker, and trial signals.

## Layout
- **`agent/`** — the deployable oncology application (FastAPI + UI, its own `docker-compose.yml`).
  This is the runnable service.
- **`docs/`** — architecture, deployment, and demo guides for the engine.

## Run
See `agent/README.md` for the full service bring-up. In brief:
```bash
cd agent
docker compose up -d      # brings up the oncology service on :8526
```

## Roadmap (this engine)
Deeper evidence integration and registration of MTB / ranking / matching as governed capabilities
in the platform registry.
