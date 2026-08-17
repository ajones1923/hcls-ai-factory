# Clinical Imaging — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/clinical-imaging/agent` · 221 Python files · 69,234 LOC · 36 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `imaging-intelligence-agent` | engine | **live** | localhost:8523 |

**UI :8523 · API :8524** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/knowledge.py`

`get_pathology_context`, `get_modality_context`, `get_anatomy_context`, `get_nim_recommendation`, `resolve_comparison_entity`, `get_comparison_context`

- **`get_pathology_context`** — Return formatted knowledge context for an imaging pathology.
- **`get_modality_context`** — Return formatted knowledge context for an imaging modality.
- **`get_anatomy_context`** — Return formatted knowledge context for an anatomical region.
- **`get_nim_recommendation`** — Return the recommended NIM workflow name for a pathology.

### `scripts/prepare_demo_data.py`

`generate_ct_chest_dicom_study`, `generate_cxr_dicom`, `precompute_workflow_results`, `precompute_radiomics`, `generate_sample_reports`, `precompute_report_nlp`

- **`generate_ct_chest_dicom_study`** — Generate a synthetic multi-slice CT chest DICOM study.
- **`generate_cxr_dicom`** — Generate a CXR DICOM from sample_cxr.dcm or sample PNG.
- **`precompute_workflow_results`** — Run all 9 workflows in mock mode, save results as JSON.
- **`precompute_radiomics`** — Generate radiomics features in mock mode.

### `src/collections.py`

`ImagingCollectionManager`

- **`ImagingCollectionManager`** — Manages 12 Imaging Milvus collections (11 owned + 1 read-only genomic).

### `src/report_parser.py`

`RadiologyReportParser`

- **`RadiologyReportParser`** — NLP parser for free-text radiology reports.

### `scripts/validate_real_data.py`

`load_metadata`, `extract_ground_truth`, `run_cxr_workflow`, `compute_metrics`, `run_validation_chest_xray`, `run_validation_pneumonia`

- **`load_metadata`** — Load metadata.json produced by download_real_data.py.
- **`extract_ground_truth`** — Convert ChestMNIST label vector to our 5-class ground truth.
- **`run_cxr_workflow`** — Run CXR classification on a single image.
- **`compute_metrics`** — Compute per-class and aggregate multi-label classification metrics.


## Dependencies

`SimpleITK>=2.3.0`, `anthropic>=0.18.0`, `apscheduler>=3.10.0`, `biopython>=1.83`, `fastapi>=0.109.0`, `highdicom>=0.22.0`, `imageio-ffmpeg>=0.4.9`, `imageio>=2.28.0`, `loguru>=0.7.0`, `lxml>=5.0.0`, `matplotlib>=3.7.0`, `monai-deploy-app-sdk>=0.6.0`, `monai>=1.3.0`, `nemoguardrails>=0.10.0`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py clinical-imaging
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

CAD-RADS is a reporting standard, not a diagnosis. The report supports a clinician's read; it does not replace it.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/clinical-imaging/agent`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
