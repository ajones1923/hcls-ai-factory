# Imaging Intelligence Agent — Demo Guide

*HCLS AI Factory — NVIDIA DGX Spark — Full Runbook*

---

## About This Guide

This guide provides detailed step-by-step procedures for every demo that can be run with the Imaging Intelligence Agent on an NVIDIA DGX Spark. It covers 29 demos organized in 7 parts, from building the system from scratch through advanced cross-modal integration.

**Each demo includes:** Purpose, Prerequisites, Step-by-step commands, Expected output, Verification checklist, and Troubleshooting tips.

**Hardware target:** NVIDIA DGX Spark — GB10 Grace Blackwell, 128 GB unified LPDDR5x, ARM64 (aarch64), NVMe storage, $4,699.

---

## Quick Reference

### Port Allocation

| Service | Host Port | Container Port | Protocol |
|---|---|---|---|
| Orthanc (DIMSE) | 4242 | 4242 | DIMSE |
| Orthanc (DICOMweb) | 8042 | 8042 | HTTP REST |
| PostgreSQL | 5432 | 5432 | TCP |
| NIM LLM | 8520 | 8000 | HTTP (OpenAI-compatible) |
| Embedding Service | 8521 | 8000 | HTTP REST |
| DICOM Listener | 8522 | 8000 | HTTP webhook |
| FHIR Publisher | 8523 | 8000 | HTTP REST |
| Agent API | 8524 | 8000 | HTTP REST |
| Streamlit Portal | 8525 | 8501 | HTTP |
| Prometheus | 9099 | 9090 | HTTP |
| Grafana | 3000 | 3000 | HTTP |
| DCGM Exporter | 9400 | 9400 | HTTP |

### Container Names

| Container | Image |
|---|---|
| imaging-orthanc | orthancteam/orthanc:24.1.2 |
| imaging-postgres | pgvector/pgvector:pg16 |
| imaging-nim-llm | nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark |
| imaging-embedding | Custom (./services/embedding_service) |
| imaging-dicom-listener | Custom (./services/dicom_listener) |
| imaging-fhir-publisher | Custom (./services/fhir_publisher) |
| imaging-agent | Custom (./agent) |
| imaging-portal | Custom (./services/portal) |
| imaging-dcgm | nvcr.io/nvidia/k8s/dcgm-exporter:3.3.5-3.4.0-ubuntu22.04 |
| imaging-prometheus | prom/prometheus:v2.48.0 |
| imaging-grafana | grafana/grafana:10.2.2 |

### Service Health Check URLs

| Service | Health Check URL |
|---|---|
| Orthanc | http://localhost:8042/system |
| NIM LLM | http://localhost:8520/v1/health/ready |
| Embedding Service | http://localhost:8521/health |
| DICOM Listener | http://localhost:8522/health |
| FHIR Publisher | http://localhost:8523/health |
| Agent API | http://localhost:8524/health |
| Streamlit Portal | http://localhost:8525/healthz |
| Prometheus | http://localhost:9099/-/healthy |
| Grafana | http://localhost:3000/api/health |

---

## Part 1: Environment Setup (Build from Scratch)

---

### Demo 0: Prerequisites & System Verification

**Purpose:** Confirm the DGX Spark hardware, software dependencies, and credentials are correctly configured before building any services.

**Prerequisites:** Physical access to DGX Spark with DGX OS installed.

#### Step 1 — Verify ARM64 Architecture

```bash
uname -m
```

**Expected output:**

```
aarch64
```

If you see `x86_64`, you are not on a DGX Spark. All containers in this project require ARM64.

#### Step 2 — Verify GPU and Memory

```bash
nvidia-smi
```

**Expected output (key lines):**

```
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 560.x.xx    Driver Version: 560.x.xx    CUDA Version: 12.x                  |
|  GPU   Name        Persistence-M | Bus-Id        Disp.A | Volatile Uncorr. ECC          |
|  0     NVIDIA GB10  On           | 00000000:01:00.0 Off |                    0          |
|                                  |                      |                                |
+-----------------------------------------------------------------------------------------+
| GPU Memory: 128 GB Unified LPDDR5x                                                      |
+-----------------------------------------------------------------------------------------+
```

#### Step 3 — Verify Docker and NVIDIA Container Toolkit

```bash
docker --version
nvidia-ctk --version
docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
```

**Expected output:** Docker version 24+, nvidia-ctk version 1.14+, and `nvidia-smi` running inside the container.

#### Step 4 — Verify NGC API Key

```bash
echo ${NGC_API_KEY:+Key is set}
```

**Expected output:**

```
Key is set
```

If blank, set it:

```bash
export NGC_API_KEY="your-ngc-api-key-here"
```

Obtain a key from https://ngc.nvidia.com/ under your account settings.

#### Step 5 — Clone Repository and Inspect Layout

```bash
git clone https://github.com/your-org/hls-imaging-agent.git
cd hls-imaging-agent
ls -la
```

**Expected output:** Repository root with `docker-compose.yml`, `main.nf`, `nextflow.config`, `maps/`, `agent/`, `services/`, `db/`, `config/`, `scripts/`, `tests/`.

#### Step 6 — Create Environment File

```bash
cp .env.example .env
```

Edit `.env` and set at minimum:

```bash
NGC_API_KEY=your-ngc-api-key
POSTGRES_USER=imaging
POSTGRES_PASSWORD=imaging_secret
POSTGRES_DB=imaging_agent
GRAFANA_USER=admin
GRAFANA_PASSWORD=changeme
```

#### Verification Checklist

- [ ] `uname -m` returns `aarch64`
- [ ] `nvidia-smi` shows GB10 with 128 GB unified memory
- [ ] Docker 24+ installed
- [ ] NVIDIA Container Toolkit installed
- [ ] NGC API key is set
- [ ] Repository cloned and `.env` file created

#### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `uname -m` shows `x86_64` | Not on DGX Spark | SSH to the correct machine |
| `nvidia-smi` not found | Driver not installed | Install DGX OS or NVIDIA driver |
| Docker permission denied | User not in docker group | `sudo usermod -aG docker $USER` then re-login |
| NGC_API_KEY empty | Not configured | Generate key at ngc.nvidia.com |

---

### Demo 1: Build & Launch All Services

**Purpose:** Build all container images and start the full 11-service stack.

**Prerequisites:** Demo 0 completed. `.env` file configured.

#### Step 1 — Build Custom Containers

```bash
docker compose build
```

This builds 6 custom services: `embedding-service`, `dicom-listener`, `fhir-publisher`, `agent`, `portal`, and the 4 MAP containers. The 5 pre-built images (orthanc, postgres, nim-llm, dcgm-exporter, prometheus, grafana) are pulled automatically.

**Expected output (final lines):**

```
[+] Building 6/6
 ✔ embedding-service  Built
 ✔ dicom-listener     Built
 ✔ fhir-publisher     Built
 ✔ agent              Built
 ✔ portal             Built
```

#### Step 2 — Launch All Services

```bash
docker compose up -d
```

**Expected output:**

```
[+] Running 11/11
 ✔ Container imaging-postgres        Healthy
 ✔ Container imaging-orthanc         Healthy
 ✔ Container imaging-nim-llm         Healthy
 ✔ Container imaging-embedding       Healthy
 ✔ Container imaging-dicom-listener  Healthy
 ✔ Container imaging-fhir-publisher  Healthy
 ✔ Container imaging-agent           Healthy
 ✔ Container imaging-portal          Started
 ✔ Container imaging-dcgm            Started
 ✔ Container imaging-prometheus      Started
 ✔ Container imaging-grafana         Started
```

#### Step 3 — Verify All Containers Healthy

```bash
docker compose ps
```

**Expected output:**

```
NAME                     IMAGE                                                    STATUS
imaging-orthanc          orthancteam/orthanc:24.1.2                              Up (healthy)
imaging-postgres         pgvector/pgvector:pg16                                  Up (healthy)
imaging-nim-llm          nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx…  Up (healthy)
imaging-embedding        imaging-embedding:latest                                Up (healthy)
imaging-dicom-listener   imaging-dicom-listener:latest                           Up (healthy)
imaging-fhir-publisher   imaging-fhir-publisher:latest                           Up (healthy)
imaging-agent            imaging-agent:latest                                    Up (healthy)
imaging-portal           imaging-portal:latest                                   Up (healthy)
imaging-dcgm             nvcr.io/nvidia/k8s/dcgm-exporter:3.3.5-3.4.0-ubunt…   Up
imaging-prometheus       prom/prometheus:v2.48.0                                 Up
imaging-grafana          grafana/grafana:10.2.2                                  Up
```

#### Step 4 — Check Each Health Endpoint

```bash
curl -s http://localhost:8042/system | python3 -m json.tool
curl -s http://localhost:8520/v1/health/ready
curl -s http://localhost:8521/health
curl -s http://localhost:8522/health
curl -s http://localhost:8523/health
curl -s http://localhost:8524/health
curl -s http://localhost:8525/healthz
curl -s http://localhost:9099/-/healthy
curl -s http://localhost:3000/api/health
```

Each should return a 200 status or JSON with `"status": "ok"`.

#### Verification Checklist

- [ ] All 11 containers are running
- [ ] Orthanc, PostgreSQL, and NIM LLM report healthy
- [ ] All custom services (embedding, dicom-listener, fhir-publisher, agent, portal) report healthy
- [ ] No port conflicts (check with `ss -tlnp | grep -E '4242|8042|5432|8520|8521|8522|8523|8524|8525|9099|3000|9400'`)

#### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| NIM LLM stays "starting" | First pull downloads model weights (~16 GB) | Wait 5-10 minutes; check `docker logs imaging-nim-llm` |
| Port already in use | Another service on that port | `ss -tlnp \| grep <port>` to find conflict; stop the other service |
| GPU access denied | Container runtime not configured | Verify `nvidia-ctk runtime configure --runtime=docker` was run |
| OOM on startup | Too many GPU services | Start services in stages: data layer first, then execution layer |

---

### Demo 2: Database Initialization

**Purpose:** Verify PostgreSQL + pgvector schema is correctly initialized with all tables, indexes, and views.

**Prerequisites:** Demo 1 completed. imaging-postgres container healthy.

#### Step 1 — Connect to PostgreSQL

```bash
docker exec -it imaging-postgres psql -U imaging -d imaging_agent
```

#### Step 2 — List Tables

```sql
\dt
```

**Expected output:**

```
              List of relations
 Schema |      Name        | Type  | Owner
--------+------------------+-------+---------
 public | embeddings       | table | imaging
 public | findings         | table | imaging
 public | measurements     | table | imaging
 public | provenance       | table | imaging
 public | series           | table | imaging
 public | studies          | table | imaging
 public | worklist_entries | table | imaging
(7 rows)
```

#### Step 3 — Verify pgvector Extension

```sql
\dx
```

**Expected output (includes):**

```
 vector | 0.7.0 | public | vector data type and ivfflat and hnsw access methods
```

#### Step 4 — Verify HNSW Index

```sql
SELECT indexname, indexdef FROM pg_indexes WHERE tablename = 'embeddings';
```

**Expected output (includes):**

```
 idx_embeddings_hnsw | CREATE INDEX idx_embeddings_hnsw ON public.embeddings
                       USING hnsw (embedding vector_cosine_ops) WITH (m='16', ef_construction='64')
```

#### Step 5 — Verify Views

```sql
\dv
```

**Expected output:**

```
          List of relations
 Schema |      Name       | Type |  Owner
--------+-----------------+------+---------
 public | active_worklist | view | imaging
 public | study_summary   | view | imaging
(2 rows)
```

#### Step 6 — Verify Empty State

```sql
SELECT COUNT(*) FROM studies;
SELECT COUNT(*) FROM findings;
SELECT COUNT(*) FROM embeddings;
```

All should return `0`.

#### Step 7 — Exit

```sql
\q
```

#### Verification Checklist

- [ ] 7 tables exist (studies, series, findings, measurements, embeddings, provenance, worklist_entries)
- [ ] pgvector extension is loaded
- [ ] HNSW index exists on embeddings table
- [ ] 2 views exist (active_worklist, study_summary)
- [ ] All tables are empty (fresh state)

---

### Demo 3: Model Downloads

**Purpose:** Download MONAI Model Zoo weights and verify NIM LLM and embedding service readiness.

**Prerequisites:** Demo 1 completed. All services running.

#### Step 1 — Download MONAI Model Zoo Weights

```bash
bash scripts/download_models.sh
```

**Expected output:**

```
Downloading DenseNet-121 (classification)...  ✓
Downloading 3D U-Net (segmentation)...        ✓
Downloading RetinaNet (detection)...           ✓
Downloading SegResNet (segmentation)...        ✓

All models downloaded to ./models/
```

#### Step 2 — Verify Model Files

```bash
ls -lh models/
```

**Expected output:**

```
-rw-r--r-- 1 user user  28M  densenet121_classification.pt
-rw-r--r-- 1 user user 120M  unet3d_segmentation.pt
-rw-r--r-- 1 user user  85M  retinanet_detection.pt
-rw-r--r-- 1 user user  95M  segresnet_segmentation.pt
```

#### Step 3 — Verify NIM LLM Readiness

```bash
curl -s http://localhost:8520/v1/models | python3 -m json.tool
```

**Expected output:**

```json
{
    "object": "list",
    "data": [
        {
            "id": "meta-llama3-8b-instruct",
            "object": "model",
            "owned_by": "nvidia"
        }
    ]
}
```

#### Step 4 — Test NIM LLM Inference

```bash
curl -s http://localhost:8520/v1/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "model": "meta-llama3-8b-instruct",
    "messages": [{"role": "user", "content": "What is Lung-RADS?"}],
    "max_tokens": 100
  }' | python3 -m json.tool
```

**Expected output:** JSON response with a chat completion about Lung-RADS classification.

#### Step 5 — Verify Embedding Service

```bash
curl -s http://localhost:8521/health
```

**Expected output:**

```json
{"status": "ok", "model": "microsoft/BiomedCLIP-PubMedBERT_256-vit_base_patch16_224"}
```

#### Verification Checklist

- [ ] All 4 MONAI model weight files present in `./models/`
- [ ] NIM LLM lists the meta-llama3-8b-instruct model
- [ ] NIM LLM responds to chat completion requests
- [ ] Embedding service reports healthy with BiomedCLIP model loaded

#### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| download_models.sh fails | Network issue or NGC auth | Check internet connectivity and NGC_API_KEY |
| NIM LLM returns 503 | Model still loading | Wait for `docker logs imaging-nim-llm` to show "ready" |
| Embedding service unhealthy | BiomedCLIP download failed | Check `docker logs imaging-embedding` for download errors |
| Disk space error | NVMe full | Check `df -h` and clean up unused Docker images |

---

### Demo 4: Orthanc DICOM Server Verification

**Purpose:** Verify the Orthanc DICOM server is operational, DICOMweb endpoints work, and the study.complete webhook pipeline fires correctly.

**Prerequisites:** Demo 1 completed. imaging-orthanc and imaging-dicom-listener containers healthy.

#### Step 1 — Check Orthanc System Info

```bash
curl -s http://localhost:8042/system | python3 -m json.tool
```

**Expected output:**

```json
{
    "ApiVersion": 23,
    "DicomAet": "IMAGING_AGENT",
    "DicomPort": 4242,
    "Name": "ImagingAgent",
    "Version": "1.12.3"
}
```

#### Step 2 — Check DICOMweb Endpoint (Empty Studies)

```bash
curl -s http://localhost:8042/dicom-web/studies | python3 -m json.tool
```

**Expected output:**

```json
[]
```

#### Step 3 — Verify DIMSE C-ECHO

```bash
docker exec imaging-orthanc echoscu localhost 4242
```

**Expected output:**

```
Association Accepted
ECHO Response: Success
Association Released
```

Alternatively, if `dcmtk` is installed locally:

```bash
echoscu localhost 4242
```

#### Step 4 — Upload a Test DICOM File

```bash
# Generate a minimal test DICOM if you don't have one
python3 -c "
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np
import tempfile, os

ds = FileDataset('test.dcm', Dataset(), preamble=b'\x00' * 128)
ds.PatientID = 'TEST001'
ds.PatientName = 'Test^Patient'
ds.StudyInstanceUID = generate_uid()
ds.SeriesInstanceUID = generate_uid()
ds.SOPInstanceUID = generate_uid()
ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2'  # CT Image Storage
ds.Modality = 'CT'
ds.StudyDate = '20260201'
ds.BodyPartExamined = 'HEAD'
ds.Rows = 512
ds.Columns = 512
ds.BitsAllocated = 16
ds.BitsStored = 12
ds.HighBit = 11
ds.PixelRepresentation = 1
ds.SamplesPerPixel = 1
ds.PhotometricInterpretation = 'MONOCHROME2'
ds.PixelData = np.random.randint(-1000, 2000, (512, 512), dtype=np.int16).tobytes()
ds.is_little_endian = True
ds.is_implicit_VR = False
ds.file_meta = pydicom.Dataset()
ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
ds.save_as('/tmp/test_ct_head.dcm')
print('Test DICOM saved to /tmp/test_ct_head.dcm')
print(f'StudyInstanceUID: {ds.StudyInstanceUID}')
"
```

#### Step 5 — Upload via DICOMweb STOW-RS

```bash
curl -X POST http://localhost:8042/dicom-web/studies \
  -H "Content-Type: application/dicom" \
  --data-binary @/tmp/test_ct_head.dcm
```

**Expected output:**

```xml
<NativeDicomModel>
  <DicomAttribute tag="00081199" vr="SQ" keyword="ReferencedSOPSequence">
    ...
  </DicomAttribute>
</NativeDicomModel>
```

#### Step 6 — Watch for Webhook (Wait StableAge = 10 Seconds)

```bash
# In a separate terminal, start watching the dicom-listener logs BEFORE uploading
docker logs -f imaging-dicom-listener
```

After 10 seconds (Orthanc `StableAge`), Orthanc fires the Lua `OnStableStudy` webhook:

**Expected output in orthanc logs:**

```
Study stable, webhook sent: <orthanc-id>
```

**Expected output in dicom-listener logs:**

```
INFO  Received study.complete webhook: orthanc_id=<orthanc-id>
INFO  Study metadata retrieved: PatientID=TEST001, Modality=CT, BodyPart=HEAD
INFO  Study inserted into database: study_id=1
INFO  Routing to workflow: ct_head_hemorrhage
```

#### Step 7 — Verify Study in Orthanc

```bash
curl -s http://localhost:8042/dicom-web/studies | python3 -m json.tool
```

**Expected output:** JSON array with 1 study entry.

#### Verification Checklist

- [ ] Orthanc system info shows AET=IMAGING_AGENT, Port=4242
- [ ] DICOMweb studies endpoint accessible
- [ ] DIMSE C-ECHO succeeds
- [ ] DICOM file uploaded via STOW-RS
- [ ] Lua webhook fires after StableAge (10 seconds)
- [ ] DICOM listener receives webhook and inserts study into database

#### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| STOW-RS returns 400 | Invalid DICOM file | Verify with `pydicom.dcmread()` first |
| Webhook never fires | Lua script not loaded | Check `docker logs imaging-orthanc` for Lua errors; verify `./config/scripts/` is mounted |
| DICOM listener doesn't receive | Network issue | Verify both containers are on `imaging-agent-network` |
| StableAge too long | Configuration | Check `orthanc.json` has `"StableAge": 10` |

---

## Part 2: Clinical Workflow Demos

---

### Demo 5: CXR Rapid Findings (Start Here — Simplest Workflow)

**Purpose:** Run an end-to-end chest X-ray classification with GradCAM heatmap generation. This is the simplest workflow (2D, single model, < 30 seconds) and the recommended starting point.

**Prerequisites:** All services running (Demo 1). DenseNet-121 model loaded (Demo 3).

**Performance target:** < 30 seconds end-to-end. Pneumothorax sensitivity > 95%.

#### Step 1 — Generate Synthetic CXR Test Data

```bash
python3 -c "
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np

ds = FileDataset('cxr.dcm', Dataset(), preamble=b'\x00' * 128)
ds.PatientID = 'CXR001'
ds.PatientName = 'Demo^CXR'
ds.StudyInstanceUID = generate_uid()
ds.SeriesInstanceUID = generate_uid()
ds.SOPInstanceUID = generate_uid()
ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.1.1'  # Digital X-Ray
ds.Modality = 'DX'
ds.StudyDate = '20260201'
ds.BodyPartExamined = 'CHEST'
ds.ViewPosition = 'PA'
ds.Rows = 2048
ds.Columns = 2048
ds.BitsAllocated = 16
ds.BitsStored = 14
ds.HighBit = 13
ds.PixelRepresentation = 0
ds.SamplesPerPixel = 1
ds.PhotometricInterpretation = 'MONOCHROME2'
ds.PixelData = np.random.randint(0, 4096, (2048, 2048), dtype=np.uint16).tobytes()
ds.is_little_endian = True
ds.is_implicit_VR = False
ds.file_meta = pydicom.Dataset()
ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
ds.save_as('/tmp/test_cxr_pa.dcm')
print('CXR test DICOM saved to /tmp/test_cxr_pa.dcm')
print(f'StudyInstanceUID: {ds.StudyInstanceUID}')
"
```

#### Step 2 — Start Watching Logs (Separate Terminal)

```bash
docker logs -f imaging-dicom-listener
```

#### Step 3 — Upload CXR via DICOMweb STOW-RS

```bash
time curl -X POST http://localhost:8042/dicom-web/studies \
  -H "Content-Type: application/dicom" \
  --data-binary @/tmp/test_cxr_pa.dcm
```

#### Step 4 — Monitor Pipeline Execution

In the dicom-listener logs, watch for:

```
INFO  Received study.complete webhook: orthanc_id=<id>
INFO  Study metadata: PatientID=CXR001, Modality=DX, BodyPart=CHEST
INFO  Routing to workflow: cxr_findings
INFO  [cxr_findings] Starting DenseNet-121 multi-label classification
INFO  [cxr_findings] Preprocessing: Resize(224,224), Normalize, EnsureChannelFirst
INFO  [cxr_findings] Inference complete (8.2s)
INFO  [cxr_findings] Results:
INFO    pneumothorax:     0.12 (negative, threshold=0.50)
INFO    consolidation:    0.87 (POSITIVE, threshold=0.40)
INFO    pleural_effusion: 0.34 (negative, threshold=0.40)
INFO    cardiomegaly:     0.62 (POSITIVE, threshold=0.50)
INFO    fracture:         0.08 (negative, threshold=0.50)
INFO  [cxr_findings] GradCAM heatmap generated for: consolidation, cardiomegaly
INFO  [cxr_findings] DICOM SR created
INFO  [cxr_findings] Findings persisted to database
INFO  [cxr_findings] Pipeline complete (12.4s total)
```

#### Step 5 — Query Findings in Database

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT f.finding_type, f.confidence, f.severity, f.is_positive
FROM findings f
JOIN studies s ON f.study_id = s.id
WHERE s.patient_id = 'CXR001'
ORDER BY f.confidence DESC;
"
```

**Expected output:**

```
 finding_type     | confidence | severity | is_positive
------------------+------------+----------+-------------
 consolidation    |       0.87 | moderate | t
 cardiomegaly     |       0.62 | moderate | t
 pleural_effusion |       0.34 | routine  | f
 pneumothorax     |       0.12 | routine  | f
 fracture         |       0.08 | routine  | f
(5 rows)
```

#### Step 6 — Check GradCAM Heatmap in Orthanc

```bash
# Find all series for the study (GradCAM stored as Secondary Capture)
curl -s "http://localhost:8042/dicom-web/studies?PatientID=CXR001" | python3 -m json.tool
```

Look for a Secondary Capture (SC) series alongside the original DX series.

#### Step 7 — Check DICOM SR

```bash
# Find SR series
curl -s "http://localhost:8042/dicom-web/series?Modality=SR&PatientID=CXR001" | python3 -m json.tool
```

#### Step 8 — Check Worklist

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT * FROM active_worklist WHERE patient_id = 'CXR001';
"
```

If consolidation or cardiomegaly exceeded their thresholds, worklist entries appear with P3 (moderate) priority.

#### Verification Checklist

- [ ] CXR uploaded successfully via STOW-RS
- [ ] Webhook fired and pipeline started within 10 seconds
- [ ] 5 finding rows in database (one per pathology)
- [ ] GradCAM heatmap generated as Secondary Capture
- [ ] DICOM SR created with coded findings
- [ ] Worklist entries created for positive findings
- [ ] Total pipeline time < 30 seconds

#### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| No pipeline execution | Webhook didn't fire | Wait full 10 seconds (StableAge); check orthanc logs |
| Model inference fails | DenseNet-121 not loaded | Check model file exists in `./models/` |
| GPU OOM | Insufficient memory | Check `nvidia-smi`; reduce batch size in config |
| GradCAM not generated | Wrong model layer target | Verify `class_layers.relu` exists in model |

---

### Demo 6: CT Head Hemorrhage Triage

**Purpose:** Detect intracranial hemorrhage, segment the bleed, measure volume and midline shift, classify urgency, and triage to the worklist.

**Prerequisites:** All services running. DenseNet-121 and 3D U-Net models loaded.

**Performance target:** < 90 seconds end-to-end. Sensitivity > 95% for hemorrhage > 5 mL.

#### Step 1 — Generate Synthetic CT Head Series

```bash
python3 -c "
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np
import os

os.makedirs('/tmp/ct_head_series', exist_ok=True)
study_uid = generate_uid()
series_uid = generate_uid()

for i in range(64):  # 64 axial slices
    ds = FileDataset(f'ct_{i:03d}.dcm', Dataset(), preamble=b'\x00' * 128)
    ds.PatientID = 'HEAD001'
    ds.PatientName = 'Demo^CTHead'
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.SOPInstanceUID = generate_uid()
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2'  # CT Image Storage
    ds.Modality = 'CT'
    ds.StudyDate = '20260201'
    ds.BodyPartExamined = 'HEAD'
    ds.StudyDescription = 'CT HEAD W/O CONTRAST'
    ds.InstanceNumber = i + 1
    ds.SliceLocation = float(i * 2.5)
    ds.ImagePositionPatient = [0.0, 0.0, float(i * 2.5)]
    ds.PixelSpacing = [0.5, 0.5]
    ds.SliceThickness = 2.5
    ds.Rows = 512
    ds.Columns = 512
    ds.BitsAllocated = 16
    ds.BitsStored = 12
    ds.HighBit = 11
    ds.PixelRepresentation = 1
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.RescaleIntercept = -1024
    ds.RescaleSlope = 1
    pixels = np.random.randint(-100, 100, (512, 512), dtype=np.int16)
    # Add a bright blob to simulate hemorrhage (HU ~60-80)
    if 20 <= i <= 35:
        pixels[200:280, 180:260] = np.random.randint(50, 80, (80, 80), dtype=np.int16)
    ds.PixelData = pixels.tobytes()
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.file_meta = pydicom.Dataset()
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.save_as(f'/tmp/ct_head_series/ct_{i:03d}.dcm')

print(f'Generated 64 CT Head slices in /tmp/ct_head_series/')
print(f'StudyInstanceUID: {study_uid}')
"
```

#### Step 2 — Upload All Slices via STOW-RS

```bash
for f in /tmp/ct_head_series/*.dcm; do
  curl -s -X POST http://localhost:8042/dicom-web/studies \
    -H "Content-Type: application/dicom" \
    --data-binary @"$f"
done
echo "All slices uploaded."
```

#### Step 3 — Wait for Pipeline (Watch Logs)

```bash
docker logs -f imaging-dicom-listener
```

**Expected pipeline stages:**

```
INFO  [ct_head_hemorrhage] Stage 1/4: DenseNet-121 classification
INFO  [ct_head_hemorrhage]   hemorrhage_detected: True (confidence=0.94)
INFO  [ct_head_hemorrhage] Stage 2/4: 3D U-Net segmentation
INFO  [ct_head_hemorrhage]   Segmentation mask generated (64 slices)
INFO  [ct_head_hemorrhage] Stage 3/4: Volumetric measurement
INFO  [ct_head_hemorrhage]   hemorrhage_volume: 12.3 mL
INFO  [ct_head_hemorrhage]   midline_shift: 3.2 mm (rightward)
INFO  [ct_head_hemorrhage]   max_thickness: 8.7 mm
INFO  [ct_head_hemorrhage] Stage 4/4: Urgency classification
INFO  [ct_head_hemorrhage]   urgency: URGENT (volume 5-30 mL)
INFO  [ct_head_hemorrhage]   priority: P2
INFO  [ct_head_hemorrhage]   routing: Neurosurgery
INFO  [ct_head_hemorrhage] Pipeline complete (47s)
```

#### Step 4 — Query Results

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT f.finding_type, f.severity, f.confidence,
       m.measurement_type, m.value, m.unit, m.flag
FROM findings f
JOIN measurements m ON m.finding_id = f.id
WHERE f.workflow = 'ct_head_hemorrhage'
ORDER BY f.id, m.measurement_type;
"
```

**Expected output:**

```
 finding_type | severity | confidence | measurement_type | value | unit | flag
--------------+----------+------------+------------------+-------+------+----------
 hemorrhage   | urgent   |       0.94 | max_thickness    |   8.7 | mm   | normal
 hemorrhage   | urgent   |       0.94 | midline_shift    |   3.2 | mm   | normal
 hemorrhage   | urgent   |       0.94 | volume           |  12.3 | mL   | elevated
```

**Urgency classification thresholds (Brain Trauma Foundation):**

| Metric | Critical | Urgent | Routine |
|---|---|---|---|
| Volume | > 30 mL | 5-30 mL | < 5 mL |
| Thickness | > 10 mm | — | <= 10 mm |
| Midline shift | > 5 mm | — | <= 5 mm |

#### Step 5 — Check Worklist Entry

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT urgency, priority, notification, routing
FROM active_worklist
WHERE patient_id = 'HEAD001';
"
```

**Expected output:**

```
 urgency | priority | notification                                | routing
---------+----------+---------------------------------------------+--------------
 urgent  | P2       | Hemorrhage detected: 12.3 mL, shift 3.2 mm | Neurosurgery
```

#### Verification Checklist

- [ ] 64 CT slices uploaded and study.complete webhook fired
- [ ] 4-stage pipeline completed: detection → segmentation → measurement → urgency
- [ ] Hemorrhage volume, midline shift, and max thickness recorded
- [ ] Urgency correctly classified (CRITICAL >30mL, URGENT 5-30mL, ROUTINE <5mL)
- [ ] Worklist entry created with correct priority (P1=critical, P2=urgent, P4=routine)
- [ ] DICOM SR with TID 1500 Measurement Report created
- [ ] Total time < 90 seconds

#### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Classification says no hemorrhage | Synthetic data may not trigger | Use real test DICOM or adjust thresholds |
| Segmentation fails | 3D U-Net not loaded | Verify `unet3d_segmentation.pt` in `./models/` |
| Volume is 0 | Empty segmentation mask | Check preprocessing transforms (windowing, spacing) |
| Wrong urgency | Threshold mismatch | Review urgency classification logic against Brain Trauma Foundation criteria |

---

### Demo 7: CT Chest Lung Nodule Tracking

**Purpose:** Detect lung nodules, measure volumes, match to prior study, calculate volume doubling time, assign Lung-RADS categories, and trigger genomics pipeline if warranted.

**Prerequisites:** All services running. RetinaNet and SegResNet models loaded.

**Performance target:** < 5 minutes multi-stage. Detection sensitivity > 90% for nodules >= 4 mm.

#### Step 1 — Upload Prior CT Chest (6 Months Earlier)

Generate and upload a prior study for the same patient:

```bash
python3 -c "
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np, os

os.makedirs('/tmp/ct_chest_prior', exist_ok=True)
study_uid = generate_uid()
series_uid = generate_uid()

for i in range(128):
    ds = FileDataset(f'ct_{i:03d}.dcm', Dataset(), preamble=b'\x00' * 128)
    ds.PatientID = 'CHEST001'
    ds.PatientName = 'Demo^CTChest'
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.SOPInstanceUID = generate_uid()
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
    ds.Modality = 'CT'
    ds.StudyDate = '20250801'  # 6 months prior
    ds.BodyPartExamined = 'CHEST'
    ds.StudyDescription = 'CT CHEST W/ CONTRAST'
    ds.InstanceNumber = i + 1
    ds.SliceLocation = float(i * 1.5)
    ds.ImagePositionPatient = [0.0, 0.0, float(i * 1.5)]
    ds.PixelSpacing = [0.7, 0.7]
    ds.SliceThickness = 1.5
    ds.Rows = 512
    ds.Columns = 512
    ds.BitsAllocated = 16
    ds.BitsStored = 12
    ds.HighBit = 11
    ds.PixelRepresentation = 1
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.RescaleIntercept = -1024
    ds.RescaleSlope = 1
    pixels = np.random.randint(-900, -500, (512, 512), dtype=np.int16)
    if 60 <= i <= 68:
        pixels[250:262, 300:312] = np.random.randint(-100, 50, (12, 12), dtype=np.int16)
    ds.PixelData = pixels.tobytes()
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.file_meta = pydicom.Dataset()
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.save_as(f'/tmp/ct_chest_prior/ct_{i:03d}.dcm')
print(f'Prior CT Chest: 128 slices in /tmp/ct_chest_prior/')
print(f'StudyInstanceUID: {study_uid}')
"

for f in /tmp/ct_chest_prior/*.dcm; do
  curl -s -X POST http://localhost:8042/dicom-web/studies \
    -H "Content-Type: application/dicom" --data-binary @"$f"
done
echo "Prior study uploaded."
```

Wait 15 seconds for the prior study pipeline to complete.

#### Step 2 — Upload Current CT Chest

```bash
python3 -c "
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np, os

os.makedirs('/tmp/ct_chest_current', exist_ok=True)
study_uid = generate_uid()
series_uid = generate_uid()

for i in range(128):
    ds = FileDataset(f'ct_{i:03d}.dcm', Dataset(), preamble=b'\x00' * 128)
    ds.PatientID = 'CHEST001'  # Same patient as prior
    ds.PatientName = 'Demo^CTChest'
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.SOPInstanceUID = generate_uid()
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
    ds.Modality = 'CT'
    ds.StudyDate = '20260201'  # Current
    ds.BodyPartExamined = 'CHEST'
    ds.StudyDescription = 'CT CHEST W/ CONTRAST'
    ds.InstanceNumber = i + 1
    ds.SliceLocation = float(i * 1.5)
    ds.ImagePositionPatient = [0.0, 0.0, float(i * 1.5)]
    ds.PixelSpacing = [0.7, 0.7]
    ds.SliceThickness = 1.5
    ds.Rows = 512
    ds.Columns = 512
    ds.BitsAllocated = 16
    ds.BitsStored = 12
    ds.HighBit = 11
    ds.PixelRepresentation = 1
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.RescaleIntercept = -1024
    ds.RescaleSlope = 1
    pixels = np.random.randint(-900, -500, (512, 512), dtype=np.int16)
    # Nodule has grown (slightly larger region)
    if 59 <= i <= 70:
        pixels[248:266, 298:316] = np.random.randint(-100, 50, (18, 18), dtype=np.int16)
    ds.PixelData = pixels.tobytes()
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.file_meta = pydicom.Dataset()
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.save_as(f'/tmp/ct_chest_current/ct_{i:03d}.dcm')
print(f'Current CT Chest: 128 slices in /tmp/ct_chest_current/')
print(f'StudyInstanceUID: {study_uid}')
"

for f in /tmp/ct_chest_current/*.dcm; do
  curl -s -X POST http://localhost:8042/dicom-web/studies \
    -H "Content-Type: application/dicom" --data-binary @"$f"
done
echo "Current study uploaded."
```

#### Step 3 — Monitor Multi-Stage Pipeline

```bash
docker logs -f imaging-dicom-listener
```

**Expected pipeline stages:**

```
INFO  [ct_chest_nodule] Stage 1/6: RetinaNet detection
INFO  [ct_chest_nodule]   Candidates found: 2
INFO  [ct_chest_nodule] Stage 2/6: SegResNet per-nodule segmentation
INFO  [ct_chest_nodule]   Nodule 1: segmented (489 mm3)
INFO  [ct_chest_nodule]   Nodule 2: segmented (112 mm3)
INFO  [ct_chest_nodule] Stage 3/6: Volumetric measurement
INFO  [ct_chest_nodule]   Nodule 1: 489 mm3, long_axis=9.8 mm
INFO  [ct_chest_nodule]   Nodule 2: 112 mm3, long_axis=5.9 mm
INFO  [ct_chest_nodule] Stage 4/6: Prior study matching
INFO  [ct_chest_nodule]   Prior study found: 2025-08-01 (184 days ago)
INFO  [ct_chest_nodule]   Registration: SyN diffeomorphic
INFO  [ct_chest_nodule]   Nodule 1 matched to prior (prior volume: 295 mm3)
INFO  [ct_chest_nodule] Stage 5/6: Volume doubling time
INFO  [ct_chest_nodule]   VDT = (184 × ln2) / ln(489/295) = 248 days
INFO  [ct_chest_nodule] Stage 6/6: Lung-RADS classification
INFO  [ct_chest_nodule]   Nodule 1: solid, 9.8 mm → Lung-RADS 4A
INFO  [ct_chest_nodule]   VDT < 400 days → upgrade to Lung-RADS 4B
INFO  [ct_chest_nodule]   Nodule 2: solid, 5.9 mm → Lung-RADS 2
INFO  [ct_chest_nodule]   GENOMICS TRIGGER: Lung-RADS 4B detected
INFO  [ct_chest_nodule] Pipeline complete (3m 12s)
```

#### Step 4 — Query Nodule Results

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT f.finding_type, f.details->>'lung_rads' AS lung_rads,
       f.confidence, f.severity,
       m.measurement_type, m.value, m.unit,
       m.prior_value, m.delta_percent
FROM findings f
JOIN measurements m ON m.finding_id = f.id
JOIN studies s ON f.study_id = s.id
WHERE s.patient_id = 'CHEST001'
  AND s.study_date = '2026-02-01'
  AND f.workflow = 'ct_chest_nodule'
ORDER BY f.id, m.measurement_type;
"
```

**Expected output:**

```
 finding_type | lung_rads | confidence | severity | measurement_type | value  | unit | prior_value | delta_percent
--------------+-----------+------------+----------+------------------+--------+------+-------------+--------------
 nodule       | 4B        |       0.91 | urgent   | doubling_time    |    248 | days |             |
 nodule       | 4B        |       0.91 | urgent   | volume           |    489 | mm3  |         295 |         65.8
 nodule       | 2         |       0.78 | routine  | volume           |    112 | mm3  |             |
```

#### Step 5 — Verify Genomics Trigger

```bash
docker logs imaging-dicom-listener | grep "GENOMICS TRIGGER"
```

**Expected output:**

```
INFO  GENOMICS TRIGGER: Lung-RADS 4B detected for study_id=<id>, triggering Parabricks pipeline
```

**VDT formula used:** `VDT = (Δt × ln2) / ln(V2/V1) = (184 × 0.693) / ln(489/295) = 248 days`

**Lung-RADS classification (solid nodules):**

| Size | Category | Growth Upgrade (VDT < 400 days) |
|---|---|---|
| < 4 mm | 1 | — |
| 4-6 mm | 2 | → 4A |
| 6-8 mm | 3 | → 4A |
| 8-15 mm | 4A | → 4B |
| >= 15 mm | 4B | — |

#### Verification Checklist

- [ ] Prior and current studies uploaded for same patient (CHEST001)
- [ ] RetinaNet detected nodule candidates
- [ ] SegResNet segmented each nodule
- [ ] Prior study retrieved and spatial registration performed
- [ ] Volume doubling time calculated correctly
- [ ] Lung-RADS assigned with growth upgrade where applicable
- [ ] Genomics pipeline triggered for Lung-RADS 4B+
- [ ] Total time < 5 minutes

---

### Demo 8: MRI Brain MS Lesion Tracking

**Purpose:** Segment MS lesions on FLAIR, count and measure them, register to prior MRI, identify new/enlarging lesions, and classify disease activity.

**Prerequisites:** All services running. 3D U-Net (FLAIR) model loaded.

**Performance target:** < 5 minutes multi-stage.

#### Step 1 — Upload Prior MRI Brain FLAIR (12 Months Earlier)

```bash
python3 -c "
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np, os

os.makedirs('/tmp/mri_brain_prior', exist_ok=True)
study_uid = generate_uid()
series_uid = generate_uid()

for i in range(96):
    ds = FileDataset(f'mr_{i:03d}.dcm', Dataset(), preamble=b'\x00' * 128)
    ds.PatientID = 'MS001'
    ds.PatientName = 'Demo^MRIBrain'
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.SOPInstanceUID = generate_uid()
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4'  # MR Image Storage
    ds.Modality = 'MR'
    ds.StudyDate = '20250201'  # 12 months prior
    ds.BodyPartExamined = 'BRAIN'
    ds.StudyDescription = 'MRI BRAIN WITH CONTRAST'
    ds.SeriesDescription = 'FLAIR'
    ds.InstanceNumber = i + 1
    ds.SliceLocation = float(i * 1.5)
    ds.ImagePositionPatient = [0.0, 0.0, float(i * 1.5)]
    ds.PixelSpacing = [1.0, 1.0]
    ds.SliceThickness = 1.5
    ds.Rows = 256
    ds.Columns = 256
    ds.BitsAllocated = 16
    ds.BitsStored = 12
    ds.HighBit = 11
    ds.PixelRepresentation = 0
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    pixels = np.random.randint(200, 800, (256, 256), dtype=np.uint16)
    # Add bright FLAIR lesion spots
    if 30 <= i <= 70:
        for cx, cy in [(80, 100), (170, 90), (120, 180)]:
            r = np.random.randint(3, 6)
            y, x = np.ogrid[-r:r+1, -r:r+1]
            mask = x**2 + y**2 <= r**2
            pixels[cy-r:cy+r+1, cx-r:cx+r+1][mask] = np.random.randint(1500, 2000)
    ds.PixelData = pixels.tobytes()
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.file_meta = pydicom.Dataset()
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.save_as(f'/tmp/mri_brain_prior/mr_{i:03d}.dcm')
print(f'Prior MRI Brain: 96 slices in /tmp/mri_brain_prior/')
"

for f in /tmp/mri_brain_prior/*.dcm; do
  curl -s -X POST http://localhost:8042/dicom-web/studies \
    -H "Content-Type: application/dicom" --data-binary @"$f"
done
echo "Prior MRI uploaded."
```

Wait 20 seconds for the prior study pipeline to complete.

#### Step 2 — Upload Current MRI Brain FLAIR

```bash
python3 -c "
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np, os

os.makedirs('/tmp/mri_brain_current', exist_ok=True)
study_uid = generate_uid()
series_uid = generate_uid()

for i in range(96):
    ds = FileDataset(f'mr_{i:03d}.dcm', Dataset(), preamble=b'\x00' * 128)
    ds.PatientID = 'MS001'  # Same patient
    ds.PatientName = 'Demo^MRIBrain'
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.SOPInstanceUID = generate_uid()
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4'
    ds.Modality = 'MR'
    ds.StudyDate = '20260201'  # Current
    ds.BodyPartExamined = 'BRAIN'
    ds.StudyDescription = 'MRI BRAIN WITH CONTRAST'
    ds.SeriesDescription = 'FLAIR'
    ds.InstanceNumber = i + 1
    ds.SliceLocation = float(i * 1.5)
    ds.ImagePositionPatient = [0.0, 0.0, float(i * 1.5)]
    ds.PixelSpacing = [1.0, 1.0]
    ds.SliceThickness = 1.5
    ds.Rows = 256
    ds.Columns = 256
    ds.BitsAllocated = 16
    ds.BitsStored = 12
    ds.HighBit = 11
    ds.PixelRepresentation = 0
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    pixels = np.random.randint(200, 800, (256, 256), dtype=np.uint16)
    # Same lesions as prior + 2 new ones
    if 30 <= i <= 70:
        for cx, cy in [(80, 100), (170, 90), (120, 180)]:
            r = np.random.randint(3, 7)  # Some slightly larger (enlarging)
            y, x = np.ogrid[-r:r+1, -r:r+1]
            mask = x**2 + y**2 <= r**2
            pixels[cy-r:cy+r+1, cx-r:cx+r+1][mask] = np.random.randint(1500, 2000)
        # 2 new lesions
        for cx, cy in [(200, 130), (60, 170)]:
            r = 4
            y, x = np.ogrid[-r:r+1, -r:r+1]
            mask = x**2 + y**2 <= r**2
            pixels[cy-r:cy+r+1, cx-r:cx+r+1][mask] = np.random.randint(1500, 2000)
    ds.PixelData = pixels.tobytes()
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.file_meta = pydicom.Dataset()
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.save_as(f'/tmp/mri_brain_current/mr_{i:03d}.dcm')
print(f'Current MRI Brain: 96 slices in /tmp/mri_brain_current/')
"

for f in /tmp/mri_brain_current/*.dcm; do
  curl -s -X POST http://localhost:8042/dicom-web/studies \
    -H "Content-Type: application/dicom" --data-binary @"$f"
done
echo "Current MRI uploaded."
```

#### Step 3 — Monitor Pipeline

**Expected pipeline stages:**

```
INFO  [mri_ms_lesion] Stage 1/5: 3D U-Net FLAIR segmentation
INFO  [mri_ms_lesion]   Segmentation mask generated (96 slices)
INFO  [mri_ms_lesion] Stage 2/5: Connected component analysis
INFO  [mri_ms_lesion]   Total lesions found: 14
INFO  [mri_ms_lesion]   Total lesion volume: 6.1 mL
INFO  [mri_ms_lesion] Stage 3/5: Prior study registration
INFO  [mri_ms_lesion]   Prior MRI found: 2025-02-01 (365 days ago)
INFO  [mri_ms_lesion]   ANTsPy SyN registration complete
INFO  [mri_ms_lesion] Stage 4/5: Change detection
INFO  [mri_ms_lesion]   Stable lesions: 11
INFO  [mri_ms_lesion]   Enlarging lesions: 1
INFO  [mri_ms_lesion]   New lesions: 2
INFO  [mri_ms_lesion] Stage 5/5: Disease activity classification
INFO  [mri_ms_lesion]   Activity: ACTIVE (1-2 new lesions)
INFO  [mri_ms_lesion] Pipeline complete (2m 48s)
```

#### Step 4 — Query Results

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT f.finding_type, f.severity,
       f.details->>'disease_activity' AS activity,
       f.details->>'new_lesion_count' AS new_lesions,
       f.details->>'enlarging_lesion_count' AS enlarging,
       m.measurement_type, m.value, m.unit,
       m.prior_value, m.delta_percent
FROM findings f
JOIN measurements m ON m.finding_id = f.id
JOIN studies s ON f.study_id = s.id
WHERE s.patient_id = 'MS001'
  AND s.study_date = '2026-02-01'
ORDER BY m.measurement_type;
"
```

**Expected output:**

```
 finding_type | severity | activity | new_lesions | enlarging | measurement_type | value | unit | prior_value | delta_percent
--------------+----------+----------+-------------+-----------+------------------+-------+------+-------------+--------------
 ms_lesion    | moderate | active   | 2           | 1         | lesion_count     |    14 | count|          12 |          16.7
 ms_lesion    | moderate | active   | 2           | 1         | volume           |   6.1 | mL   |         4.8 |          27.1
```

**Disease activity classification:**

| New/Enlarging Lesions | Activity |
|---|---|
| 0 new, 0 enlarging | Stable |
| 1-2 new or enlarging | Active |
| 3+ new or enlarging | Highly Active |

#### Verification Checklist

- [ ] Prior and current MRI uploaded for same patient (MS001)
- [ ] 3D U-Net segmentation completed on FLAIR
- [ ] Lesion count and total volume measured
- [ ] Prior study retrieved and spatial registration performed
- [ ] New and enlarging lesions identified
- [ ] Disease activity classified correctly
- [ ] Total time < 5 minutes

---

## Part 3: Agent Reasoning Demos

---

### Demo 9: Triage Agent

**Purpose:** Show the LangGraph triage agent node classifying urgency, routing to specialists, and generating notification messages.

**Prerequisites:** At least one clinical workflow completed (Demo 5, 6, 7, or 8). Agent API running.

#### Step 1 — Trigger Triage on a Study with Findings

```bash
curl -s -X POST http://localhost:8524/api/triage \
  -H "Content-Type: application/json" \
  -d '{"study_id": 1}' | python3 -m json.tool
```

**Expected output:**

```json
{
    "study_id": 1,
    "urgency": "urgent",
    "priority": "P2",
    "routing": "Neurosurgery",
    "notification": "URGENT: Intracranial hemorrhage detected. Volume 12.3 mL, midline shift 3.2 mm. Recommend immediate neurosurgical consultation.",
    "agent_trace": {
        "nodes_visited": ["triage_node"],
        "state_transitions": 1,
        "execution_time_ms": 340
    }
}
```

#### Step 2 — Verify Worklist Update

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT * FROM active_worklist ORDER BY priority, created_at;
"
```

---

### Demo 10: Longitudinal Comparison Agent

**Purpose:** Demonstrate delta computation across time points for patients with multiple studies.

**Prerequisites:** CT Chest demo (Demo 7) completed with prior and current studies.

#### Step 1 — Query Longitudinal Comparison

```bash
curl -s -X POST http://localhost:8524/api/longitudinal \
  -H "Content-Type: application/json" \
  -d '{"patient_id": "CHEST001"}' | python3 -m json.tool
```

**Expected output:**

```json
{
    "patient_id": "CHEST001",
    "studies_compared": 2,
    "time_delta_days": 184,
    "findings": [
        {
            "finding_type": "nodule",
            "current_volume_mm3": 489,
            "prior_volume_mm3": 295,
            "delta_percent": 65.8,
            "volume_doubling_time_days": 248,
            "trend": "growing",
            "lung_rads_current": "4B",
            "lung_rads_prior": "3"
        }
    ],
    "agent_trace": {
        "nodes_visited": ["longitudinal_node"],
        "execution_time_ms": 520
    }
}
```

---

### Demo 11: Population Analysis Agent

**Purpose:** Demonstrate embedding-based similarity search and cohort retrieval using pgvector.

**Prerequisites:** Multiple studies uploaded and embedded. Embedding service running.

#### Step 1 — Search for Similar Studies

```bash
curl -s -X POST http://localhost:8524/api/population \
  -H "Content-Type: application/json" \
  -d '{"study_id": 1, "top_k": 10}' | python3 -m json.tool
```

**Expected output:**

```json
{
    "query_study_id": 1,
    "similar_studies": [
        {
            "study_id": 3,
            "patient_id": "CHEST002",
            "modality": "CT",
            "study_date": "2026-01-15",
            "cosine_distance": 0.12,
            "finding_summary": "2 nodules, Lung-RADS 3"
        },
        {
            "study_id": 5,
            "patient_id": "CHEST003",
            "modality": "CT",
            "study_date": "2025-12-01",
            "cosine_distance": 0.18,
            "finding_summary": "1 nodule, Lung-RADS 4A"
        }
    ]
}
```

#### Step 2 — Direct pgvector Query

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT s.study_instance_uid, s.patient_id, s.study_date,
       e.embedding <=> (
           SELECT embedding FROM embeddings
           WHERE study_id = 1 AND level = 'study'
       ) AS distance
FROM embeddings e
JOIN studies s ON e.study_id = s.id
WHERE e.level = 'study' AND e.study_id != 1
ORDER BY distance
LIMIT 10;
"
```

---

### Demo 12: RAG-Grounded Report Generation

**Purpose:** Run the full agent pipeline (triage → longitudinal → report) with NIM LLM evidence-grounded narrative generation.

**Prerequisites:** All workflows completed. NIM LLM running.

#### Step 1 — Generate Full Report

```bash
curl -s -X POST http://localhost:8524/api/report \
  -H "Content-Type: application/json" \
  -d '{"study_id": 1}' | python3 -m json.tool
```

**Expected output (abbreviated):**

```json
{
    "study_id": 1,
    "narrative": "FINDINGS: A 12.3 mL acute right frontoparietal hemorrhage is identified with 3.2 mm rightward midline shift and maximum thickness of 8.7 mm. Per Brain Trauma Foundation guidelines, this classifies as URGENT given the volume exceeds 5 mL. Compared to similar cases in our cohort (n=8 with comparable hemorrhage volumes), early neurosurgical consultation is associated with improved outcomes.\n\nRECOMMENDATION: Immediate neurosurgical evaluation recommended.",
    "evidence_sources": [
        "findings_db: hemorrhage volume 12.3 mL, shift 3.2 mm",
        "guideline: Brain Trauma Foundation surgical criteria",
        "cohort: 8 similar cases (cosine distance < 0.2)"
    ],
    "fhir_diagnostic_report_id": "DR-2026-0001",
    "agent_trace": {
        "nodes_visited": ["triage_node", "longitudinal_node", "population_node", "report_node"],
        "total_execution_time_ms": 2840
    }
}
```

---

### Demo 13: Conditional Routing

**Purpose:** Show how the LangGraph StateGraph routes differently based on finding severity.

**Prerequisites:** Studies with CRITICAL and ROUTINE findings in the database.

#### Step 1 — CRITICAL Finding Route

```bash
curl -s -X POST http://localhost:8524/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"study_id": 1}' | python3 -c "
import json, sys
data = json.load(sys.stdin)
print('Path:', ' → '.join(data['agent_trace']['nodes_visited']))
"
```

**Expected output:**

```
Path: triage_node → longitudinal_node → population_node → report_node
```

#### Step 2 — ROUTINE Finding Route

```bash
curl -s -X POST http://localhost:8524/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"study_id": 2}' | python3 -c "
import json, sys
data = json.load(sys.stdin)
print('Path:', ' → '.join(data['agent_trace']['nodes_visited']))
"
```

**Expected output:**

```
Path: triage_node → report_node
```

Critical/urgent findings route through longitudinal and population analysis. Routine findings skip directly to report generation.

---

## Part 4: Output & Integration Demos

---

### Demo 14: DICOM SR Structured Report

**Purpose:** Inspect the TID 1500 Measurement Report content generated by the agent.

**Prerequisites:** At least one clinical workflow completed.

#### Step 1 — Find SR in Orthanc

```bash
curl -s "http://localhost:8042/dicom-web/series?Modality=SR" | python3 -m json.tool
```

#### Step 2 — Download and Parse SR

```bash
python3 -c "
from dicomweb_client import DICOMwebClient
import highdicom as hd

client = DICOMwebClient(url='http://localhost:8042/dicom-web')
series_list = client.search_for_series(search_filters={'Modality': 'SR'})

if series_list:
    study_uid = series_list[0]['0020000D']['Value'][0]
    series_uid = series_list[0]['0020000E']['Value'][0]
    instances = client.retrieve_series(study_uid, series_uid)
    sr = instances[0]
    print(f'SR Type: {sr.SOPClassUID}')
    print(f'Content Date: {sr.ContentDate}')
    # Print content tree
    for item in sr.ContentSequence:
        concept = item.ConceptNameCodeSequence[0]
        print(f'  {concept.CodeMeaning}: ', end='')
        if hasattr(item, 'MeasuredValueSequence'):
            mv = item.MeasuredValueSequence[0]
            print(f'{mv.NumericValue} {mv.MeasurementUnitsCodeSequence[0].CodeMeaning}')
        elif hasattr(item, 'TextValue'):
            print(item.TextValue)
        else:
            print('(nested)')
else:
    print('No SR found. Run a clinical workflow first.')
"
```

**Expected output:**

```
SR Type: 1.2.840.10008.5.1.4.1.1.88.34
Content Date: 20260201
  Finding: (nested)
  Volume: 12.3 mL
  Midline Shift: 3.2 mm
  Maximum Thickness: 8.7 mm
  Urgency: urgent
```

---

### Demo 15: GSPS Graphic Overlay

**Purpose:** View GradCAM heatmap overlays on original DICOM images.

**Prerequisites:** CXR Rapid Findings demo (Demo 5) completed.

#### Step 1 — Open Orthanc Explorer

Open in browser: **http://localhost:8042/app/explorer.html**

#### Step 2 — Navigate to CXR Study

Find the CXR001 patient study. You should see:
- Original DX series (the chest X-ray)
- Secondary Capture series (GradCAM heatmaps)
- GSPS series (graphic annotation overlays)

#### Step 3 — View Overlay

Click on the GSPS series to see graphic annotation contours overlaid on the original CXR. The highlighted regions correspond to GradCAM activation areas for detected pathologies.

---

### Demo 16: FHIR DiagnosticReport

**Purpose:** Inspect the FHIR R4 DiagnosticReport output with coded clinical data.

**Prerequisites:** FHIR publisher running. At least one workflow completed.

#### Step 1 — Retrieve Latest Report

```bash
curl -s http://localhost:8523/api/reports/latest | python3 -m json.tool
```

**Expected output (abbreviated):**

```json
{
    "resourceType": "Bundle",
    "type": "transaction",
    "entry": [
        {
            "resource": {
                "resourceType": "DiagnosticReport",
                "status": "final",
                "code": {
                    "coding": [{
                        "system": "http://loinc.org",
                        "code": "18748-4",
                        "display": "Diagnostic imaging study"
                    }]
                },
                "conclusion": "URGENT: Intracranial hemorrhage detected. Volume 12.3 mL.",
                "result": [
                    {"reference": "Observation/hemorrhage-volume"},
                    {"reference": "Observation/midline-shift"}
                ]
            }
        },
        {
            "resource": {
                "resourceType": "Observation",
                "id": "hemorrhage-volume",
                "code": {
                    "coding": [{
                        "system": "http://snomed.info/sct",
                        "code": "276651009",
                        "display": "Volume of hemorrhage"
                    }]
                },
                "valueQuantity": {
                    "value": 12.3,
                    "unit": "mL",
                    "system": "http://unitsofmeasure.org",
                    "code": "mL"
                }
            }
        }
    ]
}
```

---

### Demo 17: Embedding & Vector Search

**Purpose:** Directly demonstrate pgvector embedding generation and similarity queries.

**Prerequisites:** Multiple studies uploaded. Embedding service running.

#### Step 1 — Generate Embedding for a Study

```bash
curl -s -X POST http://localhost:8521/api/embed \
  -H "Content-Type: application/json" \
  -d '{"study_id": 1}' | python3 -m json.tool
```

**Expected output:**

```json
{
    "study_id": 1,
    "level": "study",
    "embedding_dim": 384,
    "model": "microsoft/BiomedCLIP-PubMedBERT_256-vit_base_patch16_224",
    "stored": true
}
```

#### Step 2 — Query Similar Studies with EXPLAIN ANALYZE

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
EXPLAIN ANALYZE
SELECT s.study_instance_uid, s.patient_id,
       e.embedding <=> (SELECT embedding FROM embeddings WHERE study_id = 1 AND level = 'study') AS distance
FROM embeddings e
JOIN studies s ON e.study_id = s.id
WHERE e.level = 'study' AND e.study_id != 1
ORDER BY distance
LIMIT 10;
"
```

**Expected output (key line):**

```
Index Scan using idx_embeddings_hnsw on embeddings e  (cost=... rows=10)
  ...
Planning Time: 0.5 ms
Execution Time: 2.3 ms
```

The HNSW index provides sub-5ms query times even with thousands of embeddings.

---

## Part 5: Portal & Monitoring Demos

---

### Demo 18: Streamlit Portal — Worklist

**Purpose:** Interactive worklist management through the web portal.

**Prerequisites:** Portal running (port 8525). Workflows completed with findings.

#### Step 1 — Open Portal

Open in browser: **http://localhost:8525**

#### Step 2 — View Active Worklist

The default view shows the active worklist sorted by priority:
- **Red** rows: P1 (critical) — immediate attention required
- **Orange** rows: P2 (urgent) — requires prompt review
- **Yellow** rows: P3 (moderate) — schedule follow-up
- **Green** rows: P4 (routine) — standard queue

Each row shows: patient ID, modality, finding type, confidence, urgency, routing target.

#### Step 3 — Drill Into a Finding

Click any worklist entry to see:
- Full finding details (type, location, confidence, severity)
- Associated measurements (volume, shift, etc.)
- Link to view images in Orthanc

#### Step 4 — Acknowledge a Finding

Click the "Acknowledge" button on a worklist entry to mark it as reviewed.

---

### Demo 19: Streamlit Portal — Study Browser

**Purpose:** Browse all processed studies with finding summaries.

#### Step 1 — Navigate to Studies Tab

Click "Studies" in the sidebar.

#### Step 2 — Browse Studies

Shows all studies with:
- Patient ID, modality, study date
- Finding count (total, critical, urgent)
- Processing status (received / processing / completed / failed)

#### Step 3 — Click Into Study Detail

Shows all findings and measurements for that study, plus links to DICOM SR and FHIR reports.

---

### Demo 20: Streamlit Portal — Agent Activity

**Purpose:** Monitor pipeline execution and model provenance.

#### Step 1 — Navigate to Activity Tab

Click "Activity" in the sidebar.

#### Step 2 — View Pipeline Runs

Shows recent provenance records:
- Workflow name
- Model ID and version
- Processing duration (ms)
- GPU memory usage (MB)
- Status (completed / failed)

---

### Demo 21: Grafana Monitoring Dashboard

**Purpose:** GPU and pipeline performance visualization.

**Prerequisites:** Grafana running (port 3000). DCGM exporter running.

#### Step 1 — Open Grafana

Open in browser: **http://localhost:3000**

Login: `admin` / `changeme`

#### Step 2 — Navigate to Imaging Agent Dashboard

Click Dashboards → Imaging Agent.

**Dashboard panels include:**

- **GPU Utilization %** — Real-time GB10 utilization from DCGM exporter
- **GPU Memory Usage** — Used vs. total (128 GB unified)
- **GPU Temperature** — Thermal monitoring
- **Inference Latency** — Histogram per workflow (p50, p95, p99)
- **Pipeline Throughput** — Studies processed per hour
- **Queue Depth** — Pending studies awaiting processing
- **Error Rate** — Failed pipeline executions

---

### Demo 22: Alerting Demo

**Purpose:** Show Prometheus alert pipeline for operational monitoring.

#### Step 1 — Check Current Alerts

```bash
curl -s http://localhost:9099/api/v1/alerts | python3 -m json.tool
```

#### Step 2 — View Alert Rules

```bash
curl -s http://localhost:9099/api/v1/rules | python3 -m json.tool
```

**Example alert rules:**

| Alert | Condition | Severity |
|---|---|---|
| GPUMemoryHigh | GPU memory > 90% for 5 min | warning |
| InferenceFailureRate | Failure rate > 5% for 10 min | critical |
| QueueBacklog | Pending studies > 50 for 15 min | warning |

---

## Part 6: Advanced Demos

---

### Demo 23: End-to-End Patient Case

**Purpose:** Upload all 4 modalities for a single patient and show unified cross-modal findings.

**Prerequisites:** All services running. All models loaded.

#### Step 1 — Upload All 4 Modalities

Use the same patient ID (`MULTI001`) for all uploads:

1. CXR (from Demo 5 pattern, change `PatientID='MULTI001'`)
2. CT Head (from Demo 6 pattern, change `PatientID='MULTI001'`)
3. CT Chest with prior (from Demo 7 pattern, change `PatientID='MULTI001'`)
4. MRI Brain with prior (from Demo 8 pattern, change `PatientID='MULTI001'`)

#### Step 2 — Wait for All Pipelines

Monitor `docker logs -f imaging-dicom-listener` until all 4 workflows complete.

#### Step 3 — Query Unified Findings

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT s.modality, f.workflow, f.finding_type, f.severity,
       f.confidence, m.measurement_type, m.value, m.unit
FROM findings f
JOIN studies s ON f.study_id = s.id
LEFT JOIN measurements m ON m.finding_id = f.id
WHERE s.patient_id = 'MULTI001'
ORDER BY
    CASE f.severity
        WHEN 'critical' THEN 1 WHEN 'urgent' THEN 2
        WHEN 'moderate' THEN 3 ELSE 4
    END,
    s.modality;
"
```

#### Step 4 — Generate Comprehensive Report

```bash
curl -s -X POST http://localhost:8524/api/patient-report \
  -H "Content-Type: application/json" \
  -d '{"patient_id": "MULTI001"}' | python3 -m json.tool
```

---

### Demo 24: Cross-Modal Genomics Trigger

**Purpose:** Demonstrate the imaging → genomics pipeline handoff when Lung-RADS 4B+ is detected.

**Prerequisites:** CT Chest demo (Demo 7) completed with Lung-RADS 4B result.

#### Step 1 — Verify Trigger in Logs

```bash
docker logs imaging-dicom-listener | grep "GENOMICS TRIGGER"
```

**Expected output:**

```
INFO  GENOMICS TRIGGER: Lung-RADS 4B detected for study_id=<id>
INFO  Parabricks pipeline trigger payload:
INFO    {"study_id": <id>, "patient_id": "CHEST001", "reason": "Lung-RADS 4B",
INFO     "workflow": "somatic_germline", "priority": "high"}
```

#### Step 2 — Inspect Trigger Payload

```bash
curl -s http://localhost:8524/api/genomics-triggers | python3 -m json.tool
```

This shows the queue of genomics pipeline triggers awaiting Parabricks execution.

---

### Demo 25: Nextflow Pipeline Orchestration

**Purpose:** Run a workflow through Nextflow and inspect the execution DAG.

**Prerequisites:** Nextflow installed. All services running.

#### Step 1 — Run Nextflow Pipeline

```bash
nextflow run main.nf \
  -profile docker \
  --study_id 1 \
  --workflow ct_head_hemorrhage \
  -with-trace \
  -with-report
```

**Expected output:**

```
N E X T F L O W  ~  version 23.10.0
Launching `main.nf` [friendly_name] DSL2 - revision: abc123

executor >  local (4)
[ab/123456] process > preprocess    [100%] 1 of 1 ✔
[cd/789012] process > inference     [100%] 1 of 1 ✔
[ef/345678] process > postprocess   [100%] 1 of 1 ✔
[gh/901234] process > persist       [100%] 1 of 1 ✔

Completed at: 2026-02-01T12:00:47
Duration    : 47s
CPU hours   : 0.01
Succeeded   : 4
```

#### Step 2 — View Execution Trace

```bash
cat trace.txt
```

Shows per-process timing, CPU, memory, and status.

#### Step 3 — Generate DAG Visualization

```bash
nextflow run main.nf -profile docker --workflow ct_head_hemorrhage -with-dag dag.html
```

Open `dag.html` in a browser to see the pipeline directed acyclic graph.

---

### Demo 26: Performance Benchmarking

**Purpose:** Measure and compare all workflow latencies against targets.

#### Step 1 — Benchmark Each Workflow

```bash
echo "=== Performance Benchmark ==="
echo ""

echo "CXR Rapid Findings:"
time curl -s -X POST http://localhost:8524/api/run-workflow \
  -d '{"study_id": <cxr_study_id>, "workflow": "cxr_findings"}' > /dev/null
echo ""

echo "CT Head Hemorrhage:"
time curl -s -X POST http://localhost:8524/api/run-workflow \
  -d '{"study_id": <ct_head_study_id>, "workflow": "ct_head_hemorrhage"}' > /dev/null
echo ""

echo "CT Chest Lung Nodule:"
time curl -s -X POST http://localhost:8524/api/run-workflow \
  -d '{"study_id": <ct_chest_study_id>, "workflow": "ct_chest_nodule"}' > /dev/null
echo ""

echo "MRI Brain MS Lesion:"
time curl -s -X POST http://localhost:8524/api/run-workflow \
  -d '{"study_id": <mri_study_id>, "workflow": "mri_ms_lesion"}' > /dev/null
```

#### Step 2 — Compare Against Targets

| Workflow | Target | Expected Actual |
|---|---|---|
| CXR Rapid Findings | < 30 seconds | ~8-12 seconds |
| CT Head Hemorrhage | < 90 seconds | ~45-60 seconds |
| CT Chest Lung Nodule | < 5 minutes | ~3-4 minutes |
| MRI Brain MS Lesion | < 5 minutes | ~2.5-3.5 minutes |

---

### Demo 27: Provenance & Reproducibility

**Purpose:** Show the complete audit trail for a finding and demonstrate deterministic re-execution.

#### Step 1 — Query Provenance

```bash
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT workflow, model_id, model_version, model_arch,
       inference_params, duration_ms, gpu_memory_mb, status,
       created_at
FROM provenance
WHERE study_id = 1
ORDER BY created_at;
"
```

**Expected output:**

```
 workflow             | model_id               | model_version | model_arch    | inference_params                    | duration_ms | gpu_memory_mb | status
----------------------+------------------------+---------------+---------------+-------------------------------------+-------------+---------------+-----------
 ct_head_hemorrhage   | hemorrhage-detect-v2.1 | 2.1.0         | DenseNet-121  | {"precision":"fp16","seed":42}      |        8200 |          3400 | completed
 ct_head_hemorrhage   | hemorrhage-seg-v1.3    | 1.3.0         | 3D U-Net      | {"precision":"fp16","seed":42}      |       28400 |         12800 | completed
```

#### Step 2 — Deterministic Re-Run

```bash
curl -s -X POST http://localhost:8524/api/reprocess \
  -H "Content-Type: application/json" \
  -d '{"study_id": 1, "workflow": "ct_head_hemorrhage"}' | python3 -m json.tool
```

With the same seed (42) and precision (fp16), outputs should be identical to the original run.

---

## Part 7: Teardown & Cleanup

---

### Demo 28: Graceful Shutdown & Data Persistence

**Purpose:** Stop all services, verify data persists across restarts, and optionally clean up completely.

#### Step 1 — Graceful Shutdown (Preserve Data)

```bash
docker compose down
```

This stops and removes containers but preserves named volumes (orthanc_data, postgres_data, nim_cache, prometheus_data, grafana_data).

#### Step 2 — Restart and Verify Persistence

```bash
docker compose up -d

# Wait for services to become healthy
sleep 30

# Verify data persisted
docker exec imaging-postgres psql -U imaging -d imaging_agent -c "
SELECT COUNT(*) AS studies, (SELECT COUNT(*) FROM findings) AS findings
FROM studies;
"
```

**Expected output:** Non-zero counts matching what was loaded before shutdown.

#### Step 3 — Full Cleanup (Remove All Data)

```bash
docker compose down -v
```

This removes all containers AND all named volumes. All data is permanently deleted.

#### Step 4 — Verify Clean State

```bash
docker volume ls | grep imaging
```

**Expected output:** No volumes listed.

---

## Appendix A — Synthetic Test Data Generation

### Generate Test DICOM for Any Modality

```python
"""Generate synthetic DICOM test data for all 4 modalities."""

import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid
import numpy as np
import os


def generate_ct_head(output_dir="/tmp/test_ct_head", num_slices=64, patient_id="TEST001"):
    """Generate synthetic CT Head with hemorrhage-like region."""
    os.makedirs(output_dir, exist_ok=True)
    study_uid = generate_uid()
    series_uid = generate_uid()

    for i in range(num_slices):
        ds = _base_dataset(patient_id, study_uid, series_uid, i)
        ds.Modality = "CT"
        ds.BodyPartExamined = "HEAD"
        ds.Rows, ds.Columns = 512, 512
        ds.PixelSpacing = [0.5, 0.5]
        ds.SliceThickness = 2.5
        ds.RescaleIntercept, ds.RescaleSlope = -1024, 1
        pixels = np.random.randint(-100, 100, (512, 512), dtype=np.int16)
        if 20 <= i <= 35:
            pixels[200:280, 180:260] = np.random.randint(50, 80, (80, 80), dtype=np.int16)
        ds.PixelData = pixels.tobytes()
        _save(ds, f"{output_dir}/ct_{i:03d}.dcm")

    print(f"CT Head: {num_slices} slices → {output_dir}")
    return study_uid


def generate_cxr(output_path="/tmp/test_cxr.dcm", patient_id="TEST002"):
    """Generate synthetic PA chest X-ray."""
    ds = _base_dataset(patient_id, generate_uid(), generate_uid(), 0)
    ds.Modality = "DX"
    ds.BodyPartExamined = "CHEST"
    ds.ViewPosition = "PA"
    ds.Rows, ds.Columns = 2048, 2048
    ds.BitsStored, ds.HighBit = 14, 13
    ds.PixelRepresentation = 0
    ds.PixelData = np.random.randint(0, 4096, (2048, 2048), dtype=np.uint16).tobytes()
    _save(ds, output_path)
    print(f"CXR: 1 image → {output_path}")


def generate_ct_chest(output_dir="/tmp/test_ct_chest", num_slices=128, patient_id="TEST003"):
    """Generate synthetic CT Chest with nodule-like region."""
    os.makedirs(output_dir, exist_ok=True)
    study_uid = generate_uid()
    series_uid = generate_uid()

    for i in range(num_slices):
        ds = _base_dataset(patient_id, study_uid, series_uid, i)
        ds.Modality = "CT"
        ds.BodyPartExamined = "CHEST"
        ds.Rows, ds.Columns = 512, 512
        ds.PixelSpacing = [0.7, 0.7]
        ds.SliceThickness = 1.5
        ds.RescaleIntercept, ds.RescaleSlope = -1024, 1
        pixels = np.random.randint(-900, -500, (512, 512), dtype=np.int16)
        if 60 <= i <= 68:
            pixels[250:262, 300:312] = np.random.randint(-100, 50, (12, 12), dtype=np.int16)
        ds.PixelData = pixels.tobytes()
        _save(ds, f"{output_dir}/ct_{i:03d}.dcm")

    print(f"CT Chest: {num_slices} slices → {output_dir}")
    return study_uid


def generate_mri_brain(output_dir="/tmp/test_mri_brain", num_slices=96, patient_id="TEST004"):
    """Generate synthetic MRI Brain FLAIR with lesion-like spots."""
    os.makedirs(output_dir, exist_ok=True)
    study_uid = generate_uid()
    series_uid = generate_uid()

    for i in range(num_slices):
        ds = _base_dataset(patient_id, study_uid, series_uid, i)
        ds.Modality = "MR"
        ds.BodyPartExamined = "BRAIN"
        ds.SeriesDescription = "FLAIR"
        ds.Rows, ds.Columns = 256, 256
        ds.PixelSpacing = [1.0, 1.0]
        ds.SliceThickness = 1.5
        ds.PixelRepresentation = 0
        pixels = np.random.randint(200, 800, (256, 256), dtype=np.uint16)
        if 30 <= i <= 70:
            for cx, cy in [(80, 100), (170, 90), (120, 180)]:
                r = np.random.randint(3, 6)
                y, x = np.ogrid[-r:r+1, -r:r+1]
                mask = x**2 + y**2 <= r**2
                pixels[cy-r:cy+r+1, cx-r:cx+r+1][mask] = np.random.randint(1500, 2000)
        ds.PixelData = pixels.tobytes()
        _save(ds, f"{output_dir}/mr_{i:03d}.dcm")

    print(f"MRI Brain: {num_slices} slices → {output_dir}")
    return study_uid


def _base_dataset(patient_id, study_uid, series_uid, instance_num):
    ds = FileDataset("temp.dcm", Dataset(), preamble=b"\x00" * 128)
    ds.PatientID = patient_id
    ds.PatientName = f"Test^{patient_id}"
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.SOPInstanceUID = generate_uid()
    ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
    ds.StudyDate = "20260201"
    ds.InstanceNumber = instance_num + 1
    ds.SliceLocation = float(instance_num * 2.5)
    ds.ImagePositionPatient = [0.0, 0.0, float(instance_num * 2.5)]
    ds.BitsAllocated = 16
    ds.BitsStored = 12
    ds.HighBit = 11
    ds.PixelRepresentation = 1
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    return ds


def _save(ds, path):
    ds.file_meta = pydicom.Dataset()
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.save_as(path)
```

---

## Appendix B — Troubleshooting Reference

| Symptom | Cause | Fix |
|---|---|---|
| Container won't start | Port conflict | `ss -tlnp \| grep <port>` to find conflict |
| NIM LLM not ready | First model download (~16 GB) | Wait 5-10 min; check `docker logs imaging-nim-llm` |
| GPU not available | NVIDIA runtime not configured | `nvidia-ctk runtime configure --runtime=docker && systemctl restart docker` |
| Webhook not firing | StableAge not elapsed | Wait 10 seconds after last instance uploaded |
| Lua script error | Script not mounted | Verify `./config/scripts:/etc/orthanc/scripts:ro` in docker-compose |
| Database connection refused | PostgreSQL not ready | Check `docker logs imaging-postgres`; verify healthcheck passes |
| Embedding service 503 | BiomedCLIP downloading | Wait for model download; check `docker logs imaging-embedding` |
| Pipeline timeout | GPU memory pressure | Check `nvidia-smi`; reduce concurrent pipelines |
| DICOM SR missing | Workflow didn't complete | Check `docker logs imaging-dicom-listener` for errors |
| Worklist empty | No positive findings | Verify findings table has rows with `is_positive = true` |
| Vector search slow | Missing HNSW index | Run `CREATE INDEX idx_embeddings_hnsw ...` from init.sql |
| ARM64 image not found | Wrong image tag | Ensure `-dgx-spark` suffix for NIM images |

---

## Appendix C — Environment Variables Reference

```bash
# .env.example — Complete environment variable reference

# ── NGC Authentication ─────────────────────────────────────
NGC_API_KEY=your-ngc-api-key-here

# ── PostgreSQL ─────────────────────────────────────────────
POSTGRES_USER=imaging
POSTGRES_PASSWORD=imaging_secret
POSTGRES_DB=imaging_agent

# ── Grafana ────────────────────────────────────────────────
GRAFANA_USER=admin
GRAFANA_PASSWORD=changeme

# ── FHIR Server (optional — external EHR integration) ─────
FHIR_SERVER_URL=http://localhost:8080/fhir

# ── Model Configuration (optional overrides) ──────────────
# EMBEDDING_MODEL=microsoft/BiomedCLIP-PubMedBERT_256-vit_base_patch16_224
# NIM_MODEL=meta-llama3-8b-instruct

# ── Pipeline Configuration (optional overrides) ───────────
# CXR_CONFIDENCE_THRESHOLD=0.50
# HEMORRHAGE_CONFIDENCE_THRESHOLD=0.50
# NODULE_CONFIDENCE_THRESHOLD=0.40
# MS_LESION_CONFIDENCE_THRESHOLD=0.50

# ── Orthanc (optional overrides) ──────────────────────────
# ORTHANC_AET=IMAGING_AGENT
# ORTHANC_STABLE_AGE=10
```

---

*HCLS AI Factory — Open Source (Apache 2.0) — NVIDIA DGX Spark*
