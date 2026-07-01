# NVIDIA FLARE Federated Learning -- Imaging Intelligence Agent

Federated learning configuration for distributed medical imaging model training across hospital sites using [NVIDIA FLARE 2.5](https://nvidia.github.io/NVFlare/). Patient imaging data never leaves the hospital -- only model weight updates are exchanged.

## Architecture

```
                           +---------------------------+
                           |     FLARE Server          |
                           |  (FL Aggregator)          |
                           |  Port 8003 (gRPC)         |
                           |  Port 8103 (Admin)        |
                           |                           |
                           |  - FedAvg Aggregation     |
                           |  - Model Persistence      |
                           |  - Cross-site Validation  |
                           |  - TensorBoard Logging    |
                           +------+----------+---------+
                                  |          |
                    Weight Diffs  |          |  Weight Diffs
                    (encrypted)   |          |  (encrypted)
                                  |          |
              +-------------------+          +-------------------+
              |                                                  |
    +---------v-----------+                        +-------------v-------+
    |   Hospital A        |                        |   Hospital B        |
    |   FLARE Client      |                        |   FLARE Client      |
    |                     |                        |                     |
    |   - Local Training  |                        |   - Local Training  |
    |   - Local Data Only |                        |   - Local Data Only |
    |   - GPU Training    |                        |   - GPU Training    |
    |   - DP (optional)   |                        |   - DP (optional)   |
    |                     |                        |                     |
    |  +---------------+  |                        |  +---------------+  |
    |  | Patient Data  |  |                        |  | Patient Data  |  |
    |  | (NEVER shared)|  |                        |  | (NEVER shared)|  |
    |  +---------------+  |                        |  +---------------+  |
    +---------------------+                        +---------------------+

    Imaging Agent Integration:
    +----------------------------------------------------------+
    |  Imaging Intelligence Agent (port 8524/8525)             |
    |  - Consumes trained global model from FLARE server       |
    |  - Uses model for inference in clinical workflows        |
    |  - CXR classification / CT segmentation / Nodule detect  |
    +----------------------------------------------------------+
```

## Federated Learning Jobs

Three imaging tasks are configured for federated training:

| Job | Model | Task | Classes | Metric | Rounds |
|-----|-------|------|---------|--------|--------|
| `cxr_classification` | DenseNet-121 | Multi-label CXR classification | 18 | AUC | 50 |
| `ct_segmentation` | 3D SegResNet | Multi-organ CT segmentation | 14 | Dice | 100 |
| `lung_nodule_detection` | 3D RetinaNet | Lung nodule detection | 2 | FROC | 80 |

## Directory Structure

```
flare/
├── README.md                           # This file
├── job_configs/
│   ├── cxr_classification/
│   │   ├── config_fed_server.json      # Server: DenseNet-121, FedAvg, AUC selection
│   │   ├── config_fed_client.json      # Client: MONAI transforms, Adam, BCEWithLogits
│   │   └── meta.json                   # 2 clients, 50 rounds, 16GB GPU
│   ├── ct_segmentation/
│   │   ├── config_fed_server.json      # Server: SegResNet, FedAvg, Dice selection
│   │   ├── config_fed_client.json      # Client: 3D MONAI transforms, AdamW, DiceCE
│   │   └── meta.json                   # 2 clients, 100 rounds, 32GB GPU
│   └── lung_nodule_detection/
│       ├── config_fed_server.json      # Server: RetinaNet, FedAvg, FROC selection
│       ├── config_fed_client.json      # Client: 3D detection transforms, SGD, Focal
│       └── meta.json                   # 2 clients, 80 rounds, 32GB GPU
├── site_configs/
│   ├── server/
│   │   └── local_config.json           # Server: gRPC 8003, admin 8103, mTLS
│   ├── site_hospital_a/
│   │   └── local_config.json           # Hospital A: data paths, privacy policies
│   └── site_hospital_b/
│       └── local_config.json           # Hospital B: data paths, privacy policies
├── provision/
│   └── project.yml                     # FLARE provisioning: participants, certs, authz
└── docker/
    └── docker-compose.flare.yml        # Full FL stack: server + 2 clients + admin + TB
```

## Prerequisites

- **NVIDIA FLARE 2.5+**: `pip install nvflare==2.5.0`
- **NVIDIA GPU**: Each client site needs at least 1 GPU with 16GB+ VRAM (32GB for 3D tasks)
- **Docker**: Docker Engine 24+ with NVIDIA Container Toolkit
- **Python 3.10+** with PyTorch 2.1+, MONAI 1.3+, torchvision
- **NGC API Key**: For pulling `nvcr.io/nvidia/nvflare` container images
- **Network**: gRPC connectivity between server and all client sites (port 8003)

### DGX Spark Compatibility

On a DGX Spark (128GB unified memory, GB10 GPU), all three containers (server + 2 simulated clients) can run on a single node for development and testing. For production, each hospital runs its own client container on-premises.

## Setup

### 1. Provision the FL Workspace

Provisioning generates startup kits with mTLS certificates, authorization policies, and startup scripts for each participant.

```bash
cd flare/

# Generate workspace (creates startup kits under /tmp/nvflare/workspace)
nvflare provision -p provision/project.yml -w provision/workspace

# The workspace will contain:
#   provision/workspace/server/         -- Server startup kit
#   provision/workspace/hospital_a/     -- Hospital A startup kit
#   provision/workspace/hospital_b/     -- Hospital B startup kit
#   provision/workspace/admin@.../      -- Admin startup kit
```

### 2. Prepare Local Data

Each hospital site must prepare its local data in the expected format. Data paths are configured in `site_configs/site_hospital_*/local_config.json`.

**CXR Classification** (CheXpert format):
```
/data/hospital_X/cxr/
├── train.csv          # Columns: Path, No Finding, Cardiomegaly, ...
├── val.csv
└── images/
    ├── patient00001/
    │   └── study1/
    │       └── view1_frontal.jpg
    ...
```

**CT Segmentation** (MONAI decathlon format):
```
/data/hospital_X/ct_seg/
├── train_data.json    # {"training": [{"image": "...", "label": "..."}]}
├── val_data.json
├── imagesTr/
│   ├── case_00001.nii.gz
│   ...
└── labelsTr/
    ├── case_00001.nii.gz
    ...
```

**Lung Nodule Detection** (LUNA16 format):
```
/data/hospital_X/lung_nodule/
├── train_data.json    # {"training": [{"image": "...", "boxes": [...]}]}
├── val_data.json
├── images/
│   ├── scan_001.nii.gz
│   ...
└── annotations/
    └── annotations.csv   # seriesuid, coordX, coordY, coordZ, diameter_mm
```

### 3. Launch the FL Environment

**Option A: Docker Compose (recommended for development)**

```bash
# Ensure the imaging-network already exists
docker network create imaging-network 2>/dev/null || true

# Start the FL stack
cd flare/docker/
docker compose -f docker-compose.flare.yml up -d

# Verify all services are running
docker compose -f docker-compose.flare.yml ps

# View server logs
docker logs -f flare-server
```

**Option B: Native (production / real multi-site)**

On each machine, run the startup script from the provisioned workspace:

```bash
# On the server machine:
cd provision/workspace/server/startup/
./start.sh

# On Hospital A's machine:
cd provision/workspace/hospital_a/startup/
./start.sh

# On Hospital B's machine:
cd provision/workspace/hospital_b/startup/
./start.sh
```

### 4. Submit a Training Job

Use the FLARE admin CLI or API to submit a federated training job.

```bash
# Enter the admin container
docker exec -it flare-admin bash

# Start the admin console
cd /workspace/admin/startup/
./fl_admin.sh

# Inside the admin console:
> check_status server
> check_status client
> submit_job /workspace/jobs/cxr_classification
> list_jobs
> abort_job <job_id>          # if needed
> shutdown client             # graceful shutdown
> shutdown server
```

**Programmatic submission** (Python API):

```python
from nvflare.fuel.flare_api.flare_api import new_secure_session

sess = new_secure_session(
    username="admin@hcls-ai-factory",
    startup_kit_dir="provision/workspace/admin@hcls-ai-factory/startup"
)

# Submit CXR classification job
job_id = sess.submit_job("job_configs/cxr_classification")
print(f"Submitted job: {job_id}")

# Monitor progress
status = sess.monitor_job(job_id, timeout=3600)
print(f"Job completed with status: {status}")

sess.close()
```

### 5. Monitor Training

**TensorBoard** (port 6006):
```bash
# If using Docker Compose, TensorBoard is already running
open http://localhost:6006

# If running natively:
tensorboard --logdir=/tmp/nvflare/log --host=0.0.0.0 --port=6006
```

**Admin CLI**:
```
> check_status server       # Shows connected clients and active jobs
> check_status client       # Shows client-side training progress
> cat server log.txt        # View server log
```

### 6. Retrieve the Trained Global Model

After a job completes, the best global model (selected by the key metric) is stored on the server.

```bash
# From the admin container or server:
ls /tmp/nvflare/models/

# Copy the model to the imaging agent for inference:
cp /tmp/nvflare/models/best_model.pt \
   /path/to/imaging_intelligence_agent/agent/data/models/

# Or use the FLARE API:
sess.download_job_result(job_id, target_dir="./results")
```

## Privacy and Security

### Data Privacy

- **Data locality**: Patient imaging data NEVER leaves the hospital site. Only model weight differences are transmitted to the server for aggregation.
- **Privacy policies**: Each site config includes `privacy.policies` with fields like `min_dataset_size` and `blocked_data_fields` to prevent accidental data leakage.

### Differential Privacy (Optional)

Each client config includes a `privacy.differential_privacy` section that can be enabled:

```json
{
  "privacy": {
    "differential_privacy": {
      "enabled": true,
      "mechanism": "gaussian",
      "epsilon": 8.0,
      "delta": 1e-5,
      "max_grad_norm": 1.0
    }
  }
}
```

When enabled, gradient clipping and calibrated Gaussian noise are applied before model updates leave the client, providing formal (epsilon, delta)-differential privacy guarantees.

| Parameter | Description | Recommended Range |
|-----------|-------------|-------------------|
| `epsilon` | Privacy budget (lower = more private) | 1.0 -- 10.0 |
| `delta` | Probability of privacy breach | 1e-5 -- 1e-7 |
| `max_grad_norm` | Gradient clipping threshold | 0.5 -- 2.0 |

### Secure Aggregation

FLARE supports homomorphic encryption (HE) for weight aggregation. The `project.yml` includes an HE builder using the CKKS scheme (supports floating point operations):

```yaml
- path: nvflare.lighter.impl.he.HEBuilder
  args:
    poly_modulus_degree: 8192
    coeff_mod_bit_sizes: [60, 40, 40]
    scale_bits: 40
    scheme: "CKKS"
```

When HE is enabled, weight updates are encrypted at each client before transmission. The server aggregates encrypted weights without ever seeing plaintext model parameters.

### Transport Security

- **mTLS**: All gRPC communication between server and clients uses mutual TLS with certificates generated during provisioning.
- **Authorization**: Role-based access control prevents unauthorized job submission or data access. Three roles are defined: `project_admin`, `org_admin`, and `lead`.

## Integration with the Imaging Agent

The federated learning pipeline produces trained global models that are consumed by the Imaging Intelligence Agent for inference.

### Workflow

```
1. Hospital sites collect and annotate local imaging data
2. FLARE trains a global model across sites (no data sharing)
3. Best global model is exported from the FLARE server
4. Model is loaded into the Imaging Agent's inference pipeline
5. Agent uses the model in clinical workflows (CXR / CT Seg / Nodule)
```

### Loading a Federated Model into the Agent

```python
import torch
from monai.networks.nets import SegResNet

# Load the federated global model
model = SegResNet(
    spatial_dims=3, in_channels=1, out_channels=14,
    init_filters=32, blocks_down=[1, 2, 2, 4], blocks_up=[1, 1, 1]
)
state_dict = torch.load("/path/to/flare/models/best_ct_seg_model.pt")
model.load_state_dict(state_dict)
model.eval()

# Use with the agent's VISTA-3D or standalone inference
```

### Port Allocation

| Service | Port | Description |
|---------|------|-------------|
| FLARE Server (gRPC) | 8003 | FL communication |
| FLARE Admin | 8103 | Admin console |
| TensorBoard | 6006 | Training visualization |
| Imaging Agent API | 8524 | REST inference |
| Imaging Agent UI | 8525 | Streamlit chat |

## Adding a New Hospital Site

1. Add the participant to `provision/project.yml`:
   ```yaml
   - name: hospital_c
     type: client
     org: hospital-c
     listening_host: 0.0.0.0
   ```

2. Create a site config directory:
   ```bash
   cp -r site_configs/site_hospital_a site_configs/site_hospital_c
   ```

3. Edit `site_configs/site_hospital_c/local_config.json`:
   - Update `client.name` to `hospital_c`
   - Update all data paths to point to Hospital C's local data
   - Adjust `resources.gpu_ids` and `resources.memory_limit_gb`

4. Re-provision to generate new certificates:
   ```bash
   nvflare provision -p provision/project.yml -w provision/workspace
   ```

5. Distribute the startup kit to Hospital C and start their client.

6. Update `meta.json` files to increase `min_clients` if desired.

## Troubleshooting

**Client cannot connect to server**:
- Verify port 8003 is open and reachable between sites
- Check that mTLS certificates were generated correctly during provisioning
- Inspect logs: `docker logs flare-hospital-a`

**Out-of-memory during training**:
- Reduce `batch_size` in the client config
- For 3D tasks, reduce `roi_size` (e.g., from [96,96,96] to [64,64,64])
- Use `sparse_update.enabled: true` with `top_k_percent: 10` to reduce communication

**Training diverges or metrics do not improve**:
- Verify data is properly formatted at each site (run local training first)
- Lower learning rate (e.g., from 0.0002 to 0.0001)
- Increase `num_rounds` and decrease local `epochs` per round
- Check for data distribution skew across sites (class imbalance)

**Job fails immediately**:
- Check `meta.json` resource requirements match client GPU availability
- Verify the model architecture matches between server and client configs
- Check the FLARE version matches (server and clients must be the same version)

## References

- [NVIDIA FLARE Documentation](https://nvidia.github.io/NVFlare/)
- [NVIDIA FLARE GitHub](https://github.com/NVIDIA/NVFlare)
- [MONAI Federated Learning Tutorial](https://github.com/Project-MONAI/tutorials/tree/main/federated_learning)
- [FLARE Medical Imaging Examples](https://github.com/NVIDIA/NVFlare/tree/main/examples/advanced/prostate)
- [CheXpert Dataset](https://stanfordmlgroup.github.io/competitions/chexpert/)
- [TotalSegmentator](https://github.com/wasserth/TotalSegmentator)
- [LUNA16 Challenge](https://luna16.grand-challenge.org/)
