# MONAI Deploy Application Packages (MAPs)

MONAI Application Packages (MAPs) wrap the HCLS AI Factory imaging workflows
into portable, DICOM-aware containers that follow clinical deployment standards.
MAPs are the standard packaging format adopted by health systems including
Mayo Clinic for deploying AI inference models into clinical PACS workflows.

## What is a MAP?

A MONAI Application Package is a containerized inference application built
with the [MONAI Deploy App SDK](https://docs.monai.io/projects/monai-deploy-app-sdk/).
Each MAP:

1. Receives DICOM studies as input (standard radiology workflow integration)
2. Selects the appropriate series based on modality and acquisition rules
3. Runs GPU-accelerated inference using the same workflow classes as the agent
4. Generates DICOM Structured Reports (SR) with coded measurements
5. Outputs JSON results for API/UI consumption

MAPs enable the same AI workflows that run in the Imaging Intelligence Agent
to be deployed directly into hospital PACS infrastructure via DICOM routing.

## Available MAPs

| MAP | Clinical Domain | Modality | Scoring System | Target Latency | Model |
|-----|----------------|----------|----------------|-----------------|-------|
| `ct_head_hemorrhage` | Neuroradiology | CT | Brain Trauma Foundation | < 90s | SegResNet |
| `cxr_rapid_findings` | Thoracic Radiology | CR/DX | Multi-label classification | < 30s | DenseNet-121 |
| `ct_chest_lung_nodule` | Thoracic Radiology | CT | Lung-RADS v2022 | < 120s | RetinaNet |
| `ct_coronary_angiography` | Cardiac Radiology | CT | CAD-RADS | < 180s | SegResNet |
| `mri_brain_ms_lesion` | Neuroradiology | MRI | MS Activity Classification | < 120s | UNEST (Large) |
| `mri_prostate_pirads` | Abdominal Radiology | MRI | PI-RADS v2.1 | < 150s | nnU-Net |
| `breast_birads` | Breast Imaging | MG/MR | BI-RADS (5th Edition) | < 60s | ResNet-50 + FPN |
| `thyroid_tirads` | Head & Neck Radiology | US | ACR TI-RADS | < 45s | EfficientNet-B4 |
| `liver_lirads` | Abdominal Radiology | CT/MRI | LI-RADS v2018 | < 180s | SegResNet + RetinaNet |

## Directory Structure

```
maps/
├── README.md
├── common/
│   └── base_imaging_map.py      # Base MONAI Deploy Application class
├── ct_head_hemorrhage/          # Emergency ICH triage (Brain Trauma Foundation)
│   ├── app.py, __main__.py, app.yaml, Dockerfile
├── cxr_rapid_findings/          # CXR multi-label classification
│   ├── app.py, __main__.py, app.yaml, Dockerfile
├── ct_chest_lung_nodule/        # Lung nodule detection (Lung-RADS v2022)
│   ├── app.py, __main__.py, app.yaml, Dockerfile
├── ct_coronary_angiography/     # Coronary CTA analysis (CAD-RADS)
│   ├── app.py, __main__.py, app.yaml, Dockerfile
├── mri_brain_ms_lesion/         # MS lesion volumetry (Activity Classification)
│   ├── app.py, __main__.py, app.yaml, Dockerfile
├── mri_prostate_pirads/         # Prostate lesion scoring (PI-RADS v2.1)
│   ├── app.py, __main__.py, app.yaml, Dockerfile
├── breast_birads/               # Breast imaging assessment (BI-RADS 5th Ed)
│   ├── app.py, __main__.py, app.yaml, Dockerfile
├── thyroid_tirads/              # Thyroid nodule scoring (ACR TI-RADS)
│   ├── app.py, __main__.py, app.yaml, Dockerfile
└── liver_lirads/                # Liver observation classification (LI-RADS v2018)
    ├── app.py, __main__.py, app.yaml, Dockerfile
```

## Relationship to Existing Workflows

MAPs wrap the same workflow classes defined in `src/workflows/`. The MAP layer
adds DICOM I/O, series selection, NIfTI conversion, and DICOM SR generation
around the existing `preprocess -> infer -> postprocess -> WorkflowResult`
pipeline. No clinical logic is duplicated.

```
DICOM Input
    │
    ▼
┌──────────────────────────────┐
│  MAP Layer (base_imaging_map)│
│  ├─ DICOM series loading     │
│  ├─ Series selection rules   │
│  ├─ DICOM → NIfTI conversion │
│  │                           │
│  │  ┌─────────────────────┐  │
│  │  │  Existing Workflow  │  │
│  │  │  preprocess()       │  │
│  │  │  infer()            │  │
│  │  │  postprocess()      │  │
│  │  │  → WorkflowResult   │  │
│  │  └─────────────────────┘  │
│  │                           │
│  ├─ DICOM SR generation      │
│  └─ JSON output              │
└──────────────────────────────┘
    │
    ▼
DICOM SR + JSON Output
```

## Building a MAP

### Prerequisites

```bash
pip install monai-deploy-app-sdk>=0.6.0
```

### Package as MAP

```bash
# Any MAP can be packaged with the same pattern:
monai-deploy package maps/<map_name>/app.py \
  -o <map_name>_map.map \
  -t <map_name>_map:1.0.0

# Examples:
monai-deploy package maps/ct_head_hemorrhage/app.py \
  -o ct_head_hemorrhage_map.map -t ct_head_hemorrhage_map:1.0.0
monai-deploy package maps/ct_chest_lung_nodule/app.py \
  -o ct_chest_lung_nodule_map.map -t ct_chest_lung_nodule_map:1.0.0
monai-deploy package maps/liver_lirads/app.py \
  -o liver_lirads_map.map -t liver_lirads_map:1.0.0
```

### Build Docker Container Directly

```bash
# From the agent/ directory — any MAP follows the same pattern:
docker build -t <map_name>_map:1.0.0 -f maps/<map_name>/Dockerfile .

# Examples:
docker build -t ct_head_hemorrhage_map:1.0.0 -f maps/ct_head_hemorrhage/Dockerfile .
docker build -t breast_birads_map:1.0.0 -f maps/breast_birads/Dockerfile .
docker build -t thyroid_tirads_map:1.0.0 -f maps/thyroid_tirads/Dockerfile .
```

## Running a MAP

### Via MONAI Deploy Runner

```bash
monai-deploy run ct_head_hemorrhage_map.map /path/to/dicom/input /path/to/output
```

### Via Docker

```bash
docker run --gpus all \
  -v /path/to/dicom:/input \
  -v /path/to/output:/output \
  ct_head_hemorrhage_map:1.0.0 --input /input --output /output --real
```

### Standalone (development/testing)

```bash
# Mock mode (no GPU required)
python -m maps.ct_head_hemorrhage --input /path/to/dicom --output /tmp/results --mock

# Real inference (requires GPU + model weights)
python -m maps.ct_head_hemorrhage --input /path/to/dicom --output /tmp/results --real
```

## Output Format

Each MAP produces two output files:

- **`result.json`** — Full `WorkflowResult` as JSON (findings, measurements,
  severity, classification, inference time)
- **`dicom_sr.dcm`** — DICOM Structured Report (TID 1500 Measurement Report)
  with coded measurements and severity classification

The DICOM SR can be sent directly to PACS via DICOM C-STORE for radiologist
review alongside the original study.

## License

Apache 2.0 — All MAPs are open-source and freely redistributable.
