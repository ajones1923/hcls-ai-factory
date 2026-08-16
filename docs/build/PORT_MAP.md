# Canonical Port Map

**Convention (adopted 2026-08-15, Adam's decision):**

> **The capability registry advertises the UI port. The API is UI + 1.**

`lib/hcls_common/capabilities.json` is the single source of truth. `scripts/validate_registry.py`
enforces the convention *and* cross-checks `health-monitor.sh`; both are in the merge gate.

---

## Allocation

| Capability | UI | API | Notes |
|---|---:|---:|---|
| genomics-engine | 5000 | — | single-service Flask portal |
| precision-intelligence-engine | 5001 | — | single-service portal |
| therapeutic-discovery-engine | 8505 | — | single-service Streamlit |
| pharmacogenomics-intelligence-agent | 8507 | 8508 | API was 8107 |
| cart-intelligence-agent | 8521 | 8522 | |
| imaging-intelligence-agent | 8523 | 8524 | UI moved 8525→8523; 8524 is the **running** cardiac API |
| precision-oncology-agent | 8526 | 8527 | already compliant |
| precision-biomarker-agent | 8528 | 8529 | UI was listed 8533 in the supervisor |
| precision-autoimmune-agent | 8531 | 8532 | already compliant |
| neurology-intelligence-agent | 8535 | 8536 | UI moved 8529→8535 (see below) |
| clinical-trial-intelligence-agent | 8538 | 8539 | API was 8538 (the UI port) |
| single-cell-intelligence-agent | 8540 | 8541 | |
| rare-disease-diagnostic-agent | 8544 | 8545 | API was 8134 |
| tuberous-sclerosis-engine | 8560 | 8561 | disease program |
| singlecell-compute | 8573 | — | compute service, no UI |
| structural-biology-engine | 8581 | 8582 | UI moved 8579→8581 |
| cardiology-intelligence-agent | 8126 | 8127 | API was 8126 (the UI port) |

**Reserved / not allocated by the convention:** 3000 Grafana · 8080 landing · 8501 RAG chat ·
8510 discovery portal · 8534 **NV-Segment-CT NIM** · 9099 Prometheus · 9100 node-exporter ·
9400 DCGM · 19530 Milvus.

---

## Why three services moved

Applying UI+1 to the old allocation produced **two collisions and two clashes with running
processes**, because four UI ports were adjacent:

| Service | Was | Now | Reason |
|---|---|---|---|
| imaging-intelligence-agent | 8525 | **8523** | its API (8526) was oncology's UI. 8524 is what the cardiac demo API actually runs on, so the pair now matches reality. |
| neurology-intelligence-agent | 8529 | **8535** | its UI was precision-biomarker's API (8528+1). Also frees **8534**, which belongs to the NV-Segment-CT NIM. |
| structural-biology-engine | 8579 | **8581** | its API (8580) is held by a running uvicorn. |

Before this, **precision-biomarker and neurology each claimed both 8528 and 8529** — they could not
run at the same time. The registry's own drift-guard could not see it, because it compared registry
entries only; the supervisor's ports were invisible to it.

---

## What the guard now catches

`validate_registry.py` fails the build on either condition. Both were **negative-tested**, not
assumed:

1. Re-introducing the neurology/biomarker overlap →
   `port convention violated on :8529 — claimed by ['neurology-intelligence-agent(ui)', 'precision-biomarker-agent(api)']`
2. Drifting a supervised port →
   `health-monitor.sh supervises 'cardiology' on :8199, which the registry does not allocate under UI/UI+1`

---

## Adding a service

1. Pick a UI port whose **UI and UI+1 are both free** — leave a gap of at least 2 from neighbours.
2. Set `endpoint: localhost:<UI>` in `capabilities.json`.
3. Add the API to `health-monitor.sh` on UI+1 (or a compose entry).
4. Run `.venv/bin/python scripts/validate_registry.py` — it will reject an adjacent allocation.

## Files changed when the convention was adopted

`lib/hcls_common/capabilities.json` (3 endpoints) · `health-monitor.sh` (8 services) ·
`core/agents/neurology/config/settings.py` · `core/agents/neurology/tests/test_settings.py` ·
neurology README / `neuro_ui.py` / DEMO_GUIDE / INDEX / DESIGN_DOCUMENT ·
`docs/HCLS_AI_FACTORY_MINDMAP.md` (6 targeted lines) · `scripts/validate_registry.py` (the guard).

> A blanket `sed` was used on the mindmap during this work and rewrote **precision-biomarker's**
> ports with neurology's. It was caught by diffing, reverted, and redone as six line-indexed edits.
> Do not run global port substitutions on platform-wide documents.
