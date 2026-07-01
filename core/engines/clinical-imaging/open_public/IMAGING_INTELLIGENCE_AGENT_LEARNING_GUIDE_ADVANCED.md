# Imaging Intelligence Agent — Learning Guide (Advanced)

> **Purpose:** Graduate-level technical deep dive into every component of the HCLS Imaging Intelligence Agent. Assumes the reader has completed the Foundations guide and is comfortable with linear algebra, statistics, Python, SQL, and software engineering.
>
> **License:** Apache 2.0 | **Author:** Adam Jones | **Date:** February 2026

---

## Table of Contents

1. [Medical Imaging Physics and Acquisition](#chapter-1--medical-imaging-physics-and-acquisition)
2. [Deep Learning Architectures for Medical Imaging](#chapter-2--deep-learning-architectures-for-medical-imaging)
3. [DICOM Deep Dive — Data Model and Networking](#chapter-3--dicom-deep-dive)
4. [GPU Architecture and CUDA Computation](#chapter-4--gpu-architecture-and-cuda-computation)
5. [Container Orchestration and Microservice Architecture](#chapter-5--container-orchestration-and-microservice-architecture)
6. [PostgreSQL + pgvector — Schema Design and Query Optimization](#chapter-6--postgresql--pgvector)
7. [Clinical Workflow Implementation — Deep Dive](#chapter-7--clinical-workflow-implementation)
8. [LangGraph Agent Architecture — Design Patterns](#chapter-8--langgraph-agent-architecture)
9. [Embedding Models and Vector Retrieval — Theory and Implementation](#chapter-9--embedding-models-and-vector-retrieval)
10. [Clinical Interoperability — FHIR, HL7, and Integration Engines](#chapter-10--clinical-interoperability)
11. [RAG Architecture and LLM Serving](#chapter-11--rag-architecture-and-llm-serving)
12. [Observability Stack — Prometheus, Grafana, and DCGM](#chapter-12--observability-stack)
13. [HCLS AI Factory — Multi-Agent Architecture](#chapter-13--hcls-ai-factory)
14. [Regulatory, Safety, and Deployment Engineering](#chapter-14--regulatory-safety-and-deployment-engineering)
15. [Discussion Questions](#discussion-questions)
16. [References](#references)

---

## Chapter 1 — Medical Imaging Physics and Acquisition

### X-Ray Attenuation and the Beer-Lambert Law

X-ray imaging exploits differential attenuation of photons through tissue. The Beer-Lambert law governs intensity reduction:

```
I = I₀ × e^(-μx)
```

Where I₀ is the incident intensity, μ is the linear attenuation coefficient (cm⁻¹), and x is the material thickness. Different tissues have different μ values — bone attenuates heavily (high μ), air barely attenuates (μ ≈ 0), and soft tissue falls between.

### Hounsfield Units (HU) and CT Windowing

CT scanners reconstruct attenuation maps calibrated to Hounsfield units:

```
HU = 1000 × (μ_tissue - μ_water) / μ_water
```

Key reference values: water = 0 HU, air = -1000 HU, dense bone = +1000 HU, acute blood = +50 to +70 HU, white matter = +20 to +30 HU.

**Windowing** maps a range of HU values to display brightness. The agent uses specific windows for each workflow:

| Window | Width (W) | Level (L) | Use |
|---|---|---|---|
| Brain | 80 | 40 | Hemorrhage detection — optimizes blood/brain contrast |
| Lung | 1500 | -600 | Lung nodule detection — shows lung parenchyma |
| Bone | 2000 | 500 | Fracture assessment |
| Soft Tissue | 400 | 50 | General abdomen/pelvis |

In DICOM, raw pixel values are converted to HU using the rescale equation:

```python
hu = pixel_value * RescaleSlope + RescaleIntercept
```

The `RescaleSlope` and `RescaleIntercept` are stored in DICOM tags (0028,1053) and (0028,1052).

### CT Reconstruction

The CT scanner acquires projections — line integrals of attenuation along ray paths. The mathematical foundation is the **Radon transform**: the set of all line integrals through a 2D function. Reconstruction recovers the original function from its Radon transform.

**Filtered Back-Projection (FBP):** The classical reconstruction algorithm. Each projection is filtered with a ramp filter (to correct for blurring inherent in back-projection) and then smeared back across the image space. FBP is fast and analytically exact for parallel-beam geometry.

**Iterative Reconstruction (IR):** Modern scanners use iterative methods (ADMIRE, ASIR-V, IMR) that model the acquisition physics, noise statistics, and image prior information. IR produces lower-noise images at reduced radiation dose but is computationally heavier.

The agent does not perform reconstruction — it receives already-reconstructed images. But understanding reconstruction artifacts (beam hardening, partial volume, motion) helps interpret model behavior on imperfect inputs.

### MRI Physics: T1/T2 Relaxation

MRI signal originates from precessing hydrogen nuclei. After RF excitation:

**T1 relaxation (longitudinal recovery):** Protons return to thermal equilibrium along the main magnetic field (B₀). T1 is the time constant for this exponential recovery:

```
Mz(t) = M₀(1 - e^(-t/T1))
```

**T2 relaxation (transverse decay):** Phase coherence is lost due to spin-spin interactions. T2 governs the exponential decay of the transverse magnetization:

```
Mxy(t) = Mxy₀ × e^(-t/T2)
```

Different tissues have characteristic T1/T2 values. CSF: long T1, long T2 (bright on T2). White matter: shorter T1 and T2 than gray matter. MS lesions: prolonged T1 and T2 (bright on T2/FLAIR).

**FLAIR Sequence Design:** FLAIR uses an inversion recovery preparation pulse timed so that CSF signal is nulled (TI chosen such that CSF's longitudinal magnetization passes through zero at readout). This makes MS lesions — which have prolonged T2 but shorter T1 than CSF — appear bright against suppressed CSF. The agent specifically requires FLAIR sequences for MS lesion segmentation.

### Voxel Anisotropy and Resampling

CT and MRI voxels are often anisotropic — the in-plane pixel spacing (e.g., 0.5 × 0.5 mm) differs from the slice thickness (e.g., 5 mm). Neural networks generally perform better on isotropic inputs.

The agent's preprocessing pipelines resample all inputs to 1.0 × 1.0 × 1.0 mm isotropic using MONAI's `Spacingd` transform with trilinear interpolation. This standardizes geometry across different scanners and protocols but introduces interpolation artifacts in the through-plane direction when slice thickness is large.

### DICOM Pixel Data Encoding

DICOM stores pixel data in the (7FE0,0010) PixelData element. Key attributes controlling interpretation:

- **Transfer Syntax:** Determines compression (1.2.840.10008.1.2 = implicit VR little-endian uncompressed; 1.2.840.10008.1.2.4.90 = JPEG 2000 lossless)
- **Photometric Interpretation:** MONOCHROME1 (white = minimum), MONOCHROME2 (white = maximum), RGB
- **Bits Allocated / Bits Stored:** Typically 16/12 or 16/16 for CT (12-bit HU range)
- **Pixel Representation:** 0 = unsigned, 1 = signed (important for HU values below zero)

Implementation with pydicom:

```python
import pydicom
import numpy as np

ds = pydicom.dcmread("slice.dcm")
pixels = ds.pixel_array  # Returns numpy array

# Convert to Hounsfield units
hu = pixels * ds.RescaleSlope + ds.RescaleIntercept
```

### Discussion Questions — Chapter 1

1. Why does the agent use a W:80 L:40 window for hemorrhage detection rather than a standard brain window?
2. What artifacts might isotropic resampling introduce when the source CT has 5 mm slice thickness?
3. How does the FLAIR inversion time determine which tissues are nulled?

---

## Chapter 2 — Deep Learning Architectures for Medical Imaging

### Convolutional Neural Networks (CNNs)

The convolution operation slides a learnable kernel (filter) across the input. For a 2D input x and kernel w of size k×k:

```
(x * w)[i,j] = Σ_m Σ_n x[i+m, j+n] × w[m, n]
```

Key hyperparameters: kernel size (typically 3×3), stride (step size between kernel applications), padding (border handling), and number of output channels (feature maps). Each layer learns multiple kernels, producing a stack of feature maps.

### DenseNet-121: Dense Connectivity

**Reference:** Huang et al., "Densely Connected Convolutional Networks," CVPR 2017.

DenseNet-121 connects each layer to every subsequent layer within a dense block. For L layers, there are L(L+1)/2 connections instead of L. Each layer receives feature maps from all preceding layers as input.

Architecture specifics:
- **4 dense blocks** with (6, 12, 24, 16) layers respectively → 6+12+24+16 = 58 dense layers + transition layers + initial conv + classifier = 121 weight layers
- **Growth rate k = 32:** Each layer produces k new feature maps. With dense connectivity, the Lth layer receives k₀ + k(L-1) input channels (where k₀ is the initial channel count)
- **Bottleneck layers:** 1×1 convolution before each 3×3 convolution to reduce computation (BN → ReLU → 1×1 Conv → BN → ReLU → 3×3 Conv)
- **Transition layers:** Between dense blocks — 1×1 convolution with compression factor θ=0.5 (halves channels) followed by 2×2 average pooling
- **Final classification:** Global average pooling → fully connected layer

The agent uses DenseNet-121 for CXR multi-label classification (5 output units with sigmoid activation) and as the initial classification gate for CT head hemorrhage detection.

### U-Net: Encoder-Decoder with Skip Connections

**Reference:** Ronneberger et al., "U-Net: Convolutional Networks for Biomedical Image Segmentation," MICCAI 2015.

U-Net has a symmetric encoder-decoder architecture:

- **Encoder (contracting path):** Repeated application of two 3×3 convolutions (each followed by batch normalization and ReLU), then 2×2 max pooling for downsampling. Each downsampling step doubles the feature channels.
- **Decoder (expanding path):** 2×2 transposed convolution for upsampling, concatenation with the corresponding encoder feature maps (skip connections), then two 3×3 convolutions.
- **Skip connections:** Copy feature maps from encoder to decoder at each resolution level. This preserves spatial detail that the encoder's pooling operations would otherwise lose.

The 3D U-Net extension processes volumetric data (5D input: batch × channels × depth × height × width). Memory requirements scale cubically with spatial dimensions — a 512×512×300 input at full resolution would require enormous GPU memory. The agent addresses this with:
- Isotropic resampling to 1mm (reducing dimensions)
- Sliding window inference (processing overlapping patches)
- Reduced channel counts for the first layers

The agent's hemorrhage segmentation model:
```python
UNet(spatial_dims=3, in_channels=1, out_channels=2,
     channels=(16, 32, 64, 128, 256), strides=(2, 2, 2, 2),
     num_res_units=2, norm="batch")
```

### SegResNet: Residual Blocks in Segmentation

**Reference:** Myronenko, "3D MRI Brain Tumor Segmentation Using Autoencoder Regularization," BrainLes@MICCAI 2018.

SegResNet introduces residual blocks (from ResNet) into a U-Net-like architecture. Each block computes:

```
output = F(x) + x
```

Where F(x) is the block's learned transformation. The identity shortcut (+x) addresses the vanishing gradient problem — gradients flow directly through the skip connection during backpropagation, enabling training of deeper networks.

The agent uses SegResNet for per-nodule lung segmentation:
```python
SegResNet(spatial_dims=3, in_channels=1, out_channels=2,
          init_filters=16, blocks_down=(1, 2, 2, 4),
          blocks_up=(1, 1, 1), norm="batch")
```

### RetinaNet: Feature Pyramid Network + Focal Loss

**Reference:** Lin et al., "Focal Loss for Dense Object Detection," ICCV 2017.

RetinaNet addresses the class imbalance problem in detection (overwhelmingly many background regions vs. few objects). Two key innovations:

**Feature Pyramid Network (FPN):** Builds a multi-scale feature pyramid from the backbone network (ResNet-50). The pyramid has levels P3 through P7, with each level operating at a different spatial resolution. Objects of different sizes are detected at different pyramid levels.

**Focal Loss:** Modifies the standard cross-entropy loss to down-weight easy examples:

```
FL(p_t) = -α_t (1 - p_t)^γ × log(p_t)
```

Where p_t is the predicted probability for the correct class, γ is the focusing parameter (typically γ=2), and α is a class-balancing weight. When a sample is correctly classified with high confidence (p_t → 1), the modulating factor (1-p_t)^γ → 0, reducing the loss contribution. This focuses training on hard, misclassified examples.

**Anchor design:** At each FPN level, anchors of different scales and aspect ratios tile the feature map. The agent uses sizes (4, 6, 8, 12, 16) and aspect ratios (0.5, 1.0, 2.0) for nodule detection. Non-maximum suppression (NMS) with threshold 0.3 removes duplicate detections.

### SwinUNETR: Transformer-Based Segmentation

**Reference:** Hatamizadeh et al., "Swin UNETR: Swin Transformers for Semantic Segmentation of Brain Tumors in MRI Images," BrainLes@MICCAI 2021.

SwinUNETR uses a Swin Transformer encoder (shifted window self-attention) with a CNN decoder:

- **Patch embedding:** Input volume is partitioned into non-overlapping 3D patches (e.g., 2×2×2 voxels → tokens)
- **Shifted window attention:** Self-attention computed within local windows, with alternate layers shifting the window partition to enable cross-window connections
- **Hierarchical features:** 4 stages with downsampling via patch merging — analogous to the CNN encoder's pooling
- **Decoder:** CNN-based with skip connections from each transformer stage

SwinUNETR is available in the MONAI Model Zoo and can be used as a drop-in replacement for the U-Net in segmentation workflows when transformer-based spatial modeling is beneficial.

### Loss Functions

| Loss | Formula | When to Use |
|---|---|---|
| Cross-Entropy | -Σ y_i log(p_i) | Classification, per-pixel classification |
| Dice Loss | 1 - 2\|A∩B\|/(\|A\|+\|B\|) | Segmentation — handles class imbalance |
| Focal Loss | -α(1-p_t)^γ log(p_t) | Detection — focuses on hard examples |
| Combined (Dice+CE) | λ₁×Dice + λ₂×CE | Segmentation — benefits from both |

The agent uses Dice + cross-entropy combined loss for segmentation training and focal loss for RetinaNet detection.

### Data Augmentation for Medical Imaging

Medical imaging datasets are small relative to natural image datasets. Standard augmentations:

- **Affine:** Random rotation (±15°), scaling (0.85-1.15), translation (±10 pixels)
- **Elastic deformation:** Random B-spline displacement field — simulates anatomical variability
- **Intensity augmentation:** Random brightness/contrast shifts, noise injection (Gaussian, Rician for MRI)
- **Flipping:** Left-right mirroring (valid for bilateral anatomy but must be avoided for laterality-dependent tasks)

MONAI provides all of these via the transforms pipeline (e.g., `RandAffined`, `RandElasticd`, `RandGaussianNoised`).

### Transfer Learning

Pre-trained weights accelerate convergence and improve performance on small datasets. Two strategies:

1. **ImageNet pre-training:** General visual features. Wide availability. But significant domain gap from medical imaging (color photos vs. grayscale, vastly different subjects).
2. **Medical pre-training:** Weights trained on large medical image datasets (e.g., CheXpert for CXR, BTCV for CT). Closer domain match. MONAI Model Zoo provides medical pre-trained weights for most architectures.

### MONAI Framework Patterns

MONAI (Medical Open Network for AI, Apache 2.0) provides:

- **Dictionary-based transforms:** Operations keyed by data role (`keys=["image", "label"]`) enabling matched augmentation of image-label pairs
- **Sliding window inference:** For volumes larger than GPU memory — processes overlapping patches, then stitches outputs with weighted averaging in overlap regions
- **Dataset classes:** `CacheDataset` (pre-caches transforms in memory), `PersistentDataset` (caches to disk)

### Discussion Questions — Chapter 2

1. Derive the total number of parameters in the classification head of DenseNet-121 with growth rate k=32 given 4 dense blocks of (6,12,24,16) layers.
2. Why does focal loss with γ=2 outperform standard cross-entropy for lung nodule detection?
3. What is the memory complexity of 3D sliding window inference with window size W³ and overlap ratio r?

---

## Chapter 3 — DICOM Deep Dive

### Information Object Definitions (IODs)

DICOM organizes data in a four-level hierarchy: Patient → Study → Series → Instance. Each level has specific Information Object Definitions (IODs) that define required and optional attributes.

**Composite IODs** (most common) include attributes from multiple Information Entities. A CT Image IOD includes Patient IE, Study IE, Series IE, Frame of Reference IE, Equipment IE, and Image IE attributes — all in a single DICOM file.

**Normalized IODs** represent a single entity (e.g., Modality Worklist) and are accessed via N-services (N-GET, N-SET, N-CREATE) rather than C-services.

### Value Representations (VR)

Each DICOM attribute has a Value Representation — its data type:

| VR | Name | Example |
|---|---|---|
| PN | Person Name | DOE^JOHN^A^^DR |
| DA | Date | 20260115 |
| TM | Time | 143205.123456 |
| UI | Unique Identifier | 1.2.840.113619.2.55.1 |
| DS | Decimal String | "25.3" (stored as text) |
| IS | Integer String | "300" |
| SQ | Sequence | Nested dataset(s) |
| OW/OB | Other Word/Byte | Pixel data |

Transfer syntaxes determine byte ordering and whether VRs are encoded explicitly or implicitly. DICOM's default (1.2.840.10008.1.2) is Implicit VR Little Endian, but Explicit VR Little Endian (1.2.840.10008.1.2.1) is preferred because the VR is stored in the data stream, eliminating dictionary lookup requirements.

### DIMSE Services

DICOM Message Service Element (DIMSE) services operate over TCP connections:

| Service | Direction | Purpose |
|---|---|---|
| C-ECHO | SCU → SCP | Verify connectivity ("ping") |
| C-STORE | SCU → SCP | Send a DICOM instance |
| C-FIND | SCU → SCP | Query for matching instances |
| C-MOVE | SCU → SCP → 3rd party | Request SCP to send instances to a third-party AE |
| C-GET | SCU → SCP | Request SCP to send instances back to the SCU |

Before any DIMSE operation, an **Association negotiation** occurs. The SCU proposes presentation contexts (SOP Class + Transfer Syntax pairs), and the SCP accepts or rejects each. Only accepted contexts can be used.

The agent's Orthanc server accepts C-STORE for all CT, MR, CR, and DX SOP Classes, enabling scanners and PACS to push images directly.

### DICOMweb REST API

DICOMweb maps DICOM operations to HTTP REST:

**STOW-RS (Store):** POST multipart/related request with DICOM instances in the request body. Content-Type per part: `application/dicom`.

**WADO-RS (Retrieve):** GET request with study/series/instance UIDs in the URL path. Supports metadata retrieval (JSON or XML), bulk data retrieval, and rendered frame retrieval.

**QIDO-RS (Search):** GET request with query parameters mapped from DICOM tag names. Supports `limit`, `offset`, `includefield`, and standard DICOM matching (exact, wildcard, range for dates/times).

### DICOM SR Template TID 1500

The agent's findings are encoded as DICOM SR Measurement Reports following template TID 1500. The content tree structure:

```
Measurement Report (root)
├── Language of Content (EN)
├── Observation Context
│   ├── Observer Type (Device)
│   └── Device Observer (Imaging Intelligence Agent v2.1)
├── Image Library (reference to source images)
└── Measurement Group(s)
    ├── Tracking Identifier
    ├── Finding (coded concept, e.g., SCT:50960005 = Hemorrhage)
    ├── Measurement(s)
    │   ├── Concept Name (e.g., DCM:118565006 = Volume)
    │   ├── Measured Value (25.3)
    │   └── Units (UCUM: mL)
    └── Qualitative Evaluation(s)
        └── Severity (CRITICAL)
```

Each measurement includes **SCOORD** (Spatial Coordinates) references linking the measurement to specific image frames and regions.

The agent uses **highdicom** to construct SRs:
```python
import highdicom as hd

measurement = hd.sr.Measurement(
    name=hd.sr.CodedConcept("118565006", "SCT", "Volume"),
    value=25.3,
    unit=hd.sr.CodedConcept("mL", "UCUM", "milliliter"),
)
```

### Orthanc Architecture

Orthanc (GPLv3) is a lightweight DICOM server with:
- **SQLite metadata index:** Fast lookups by UID, patient ID, date
- **On-disk storage:** Each DICOM instance stored as a separate file, organized by hash
- **Lua scripting hooks:** Event-driven callbacks (OnStoredInstance, OnStableStudy) — the agent's trigger mechanism
- **REST API:** Full CRUD on patients, studies, series, instances, plus DICOMweb endpoints
- **Plugin architecture:** DICOMweb, PostgreSQL index, authorization modules

The agent's Lua callback fires when a study has been stable (no new instances received) for `StableAge` seconds (configured as 10):

```lua
function OnStableStudy(studyId, tags, metadata)
    -- Port 8000 is the container-internal port; mapped to host 8522 in docker-compose
    local url = "http://dicom-listener:8000/webhook/study-complete"
    local body = '{"orthanc_id": "' .. studyId .. '"}'
    HttpPost(url, body, {["Content-Type"] = "application/json"})
end
```

### Discussion Questions — Chapter 3

1. Why does the agent use StableAge=10 seconds instead of processing each instance immediately?
2. Describe the Association negotiation required for Orthanc to accept a JPEG 2000 compressed CT series.
3. How would you extend the DICOM SR template to include lesion-level embedding references?

---

## Chapter 4 — GPU Architecture and CUDA Computation

### GPU Hardware Architecture

A modern NVIDIA GPU consists of an array of **Streaming Multiprocessors (SMs)**. Each SM contains:
- CUDA cores (FP32 units)
- Tensor Cores (matrix multiply-accumulate units optimized for deep learning)
- Warp schedulers (dispatch units that execute threads in groups of 32 called "warps")
- Shared memory and L1 cache (configurable split)
- Register file

The thread hierarchy: a **grid** contains **blocks**, each block contains **threads**. The hardware maps blocks to SMs and threads to CUDA cores within an SM. Threads within a block can synchronize and share data via shared memory.

### Memory Hierarchy

From fastest to slowest:
1. **Registers** (~1 cycle latency, per-thread)
2. **Shared memory / L1 cache** (~20-30 cycles, per-SM, configurable 48-128 KB)
3. **L2 cache** (~200 cycles, shared across all SMs)
4. **Global memory** (HBM2e on data center GPUs, LPDDR5x on DGX Spark, ~400-600 cycles)

On the DGX Spark, the 128 GB LPDDR5x is unified — the same physical memory is accessible by both the Grace CPU and the GB10 GPU. This eliminates explicit PCIe transfers but introduces coherence overhead when both processors access the same pages.

### Grace Blackwell GB10

The DGX Spark's GB10 is a Blackwell-generation GPU with:
- Compute capability 10.0+
- 5th-generation Tensor Cores (FP4, FP8, FP16, BF16, TF32 support)
- Hardware-accelerated sparsity (2:4 structured sparsity for 2× throughput)
- Unified memory with Grace CPU via NVLink-C2C (coherent interconnect, not PCIe)

The NVLink-C2C interconnect between Grace and Blackwell provides ~900 GB/s bidirectional bandwidth — far exceeding PCIe Gen5 (128 GB/s). Combined with unified memory addressing, this means the GPU can directly access CPU memory allocations without explicit copy commands, though access patterns affect performance.

### cuDNN Convolution Algorithms

cuDNN offers multiple algorithms for convolution, each optimal for different tensor shapes:

| Algorithm | Approach | Best For |
|---|---|---|
| Implicit GEMM | Converts conv to matrix multiply | General purpose, reliable |
| FFT | Transform to frequency domain, multiply, transform back | Large kernels, large batches |
| Winograd | Minimal filtering algorithm (fewer multiplications) | 3×3 and 5×5 kernels |

cuDNN's **autotuning** (`torch.backends.cudnn.benchmark = True`) benchmarks all algorithms on the first forward pass and selects the fastest for each layer configuration. This is recommended for the agent since input shapes are fixed per workflow.

### TensorRT Optimization

TensorRT converts a trained PyTorch/ONNX model into an optimized inference engine:

1. **Graph optimization:** Operator fusion (Conv + BN + ReLU → single kernel), constant folding, dead code elimination
2. **Precision calibration:** FP32 → FP16 (2× throughput, minimal accuracy loss) or INT8 (4× throughput, requires calibration dataset to minimize quantization error)
3. **Kernel autotuning:** Selects optimal CUDA kernels per layer based on the specific GPU
4. **Memory optimization:** Layer-level scratch space allocation, concurrent kernel execution planning

For the agent's workflows, FP16 inference provides the best tradeoff — nearly 2× speedup with negligible accuracy impact on segmentation tasks. INT8 is viable for classification (DenseNet-121) but requires careful validation for segmentation models.

### GPUDirect Storage

GPUDirect Storage (GDS) eliminates the CPU from the I/O path:

**Without GDS:** NVMe → OS page cache → CPU memory → PCIe/NVLink → GPU memory (3 copies)
**With GDS:** NVMe → GPU memory via DMA (1 copy, zero CPU bounce buffers)

GDS uses the cuFile API:
```c
cuFileDriverOpen();
cuFileHandleRegister(&fh, &descr);
cuFileRead(fh, devPtr, size, file_offset, 0);
```

For the agent, GDS accelerates DICOM loading when studies are stored on local NVMe — particularly beneficial for large CT volumes (500+ MB) and batch processing scenarios.

### Memory Budget Analysis

For a given model and input, total GPU memory = model parameters + activations + workspace:

```
Model memory ≈ num_parameters × bytes_per_param (4 for FP32, 2 for FP16)
Activation memory ≈ batch_size × sum(output_size_per_layer × bytes_per_element)
```

The agent's 3D U-Net for hemorrhage segmentation (~2M parameters) at FP16 requires ~4 MB for weights. But activations for a 256×256×256 input volume can require 2-4 GB. On the DGX Spark's 128 GB unified memory, this is comfortable — but memory must be shared with the NIM LLM service (~16 GB for Llama 3 8B), the embedding model (~500 MB), and the OS.

### ARM64 Considerations

The Grace CPU is ARM64 (aarch64). Critical implications:

- **Python wheels:** Most scientific packages (numpy, scipy, PyTorch, MONAI) now provide ARM64 wheels, but some niche packages may require compilation from source
- **NEON SIMD:** ARM's SIMD instruction set (128-bit vectors). Performance-sensitive C extensions must be compiled with NEON support
- **SVE (Scalable Vector Extension):** Grace supports SVE2 with variable-width vectors. Not yet widely leveraged by Python ecosystem
- **Docker base images:** Must use ARM64 variants (e.g., `--platform linux/arm64`)

### Discussion Questions — Chapter 4

1. Calculate the theoretical peak FP16 throughput for the GB10 given its SM count and tensor core capabilities.
2. What is the expected speedup from TensorRT FP16 optimization for a 3D U-Net inference compared to PyTorch eager mode FP32?
3. Under what circumstances would GPUDirect Storage not improve performance on the DGX Spark?

---

## Chapter 5 — Container Orchestration and Microservice Architecture

### OCI Container Specification

Containers follow the Open Container Initiative (OCI) specification:

- **Image manifest:** JSON document listing layers, config, and platform
- **Config:** Entry point, environment variables, exposed ports, volumes
- **Layers:** Filesystem diffs stacked via union mount (overlay2 driver on Linux)
- **Rootfs:** The combined filesystem visible to the container

Each layer is content-addressable (SHA-256 hash). Identical layers are shared across images, reducing disk usage.

### Multi-Stage Builds

The agent's Dockerfiles use multi-stage builds to minimize image size:

```dockerfile
# Stage 1: Build dependencies
FROM python:3.11-slim AS builder
RUN pip install --no-cache-dir --target=/deps -r requirements.txt

# Stage 2: Runtime image
FROM python:3.11-slim
COPY --from=builder /deps /usr/local/lib/python3.11/site-packages
COPY app.py .
```

Build stage artifacts (compilers, headers, pip cache) are excluded from the final image.

### ARM64 Multi-Architecture Builds

For DGX Spark compatibility:

```bash
docker buildx create --use --name multiarch
docker buildx build --platform linux/arm64 -t imaging-agent:latest .
```

Alternatively, build natively on the DGX Spark itself (no cross-compilation needed). QEMU user-space emulation (`docker run --privileged multiarch/qemu-user-static --reset -p yes`) enables ARM64 builds on x86 hosts but is significantly slower.

### The Agent's Service Topology

Dependency graph (arrows indicate "depends on"):

```
orthanc ←── dicom-listener ←── agent ←── portal
                                 ↑
postgres ←── embedding-service   │
    ↑                            │
    └── fhir-publisher           │
                                 │
nim-llm ─────────────────────────┘

dcgm-exporter ←── prometheus ←── grafana
```

**Startup ordering:** docker-compose `depends_on` with `condition: service_healthy` ensures services start only after their dependencies pass health checks. The longest startup path is: postgres → nim-llm → agent → portal (NIM model loading dominates at ~60-120 seconds).

**Failure domains:** If postgres fails, all data services are affected. If nim-llm fails, report generation degrades but triage still works. If orthanc fails, no new studies can be ingested but existing data remains queryable. This separation limits blast radius.

### Health Check Patterns

Three health check approaches in the agent:

| Pattern | Implementation | Used By |
|---|---|---|
| HTTP GET | `curl -f http://localhost:PORT/health` | All custom services (8520-8525) |
| TCP Socket | pg_isready utility | PostgreSQL |
| DIMSE echo | curl to Orthanc REST `/system` endpoint | Orthanc |

Health check parameters: `interval=30s, timeout=10s, retries=3` for most services. Shorter interval (10s) for PostgreSQL since it is a critical dependency.

### Container Security

Production hardening measures:
- **Non-root execution:** `USER 1000:1000` in Dockerfile — prevents container escape privilege escalation
- **Read-only filesystem:** `read_only: true` in compose with explicit tmpfs for /tmp
- **Capability dropping:** `cap_drop: [ALL]` with selective `cap_add` (only what is needed)
- **No new privileges:** `security_opt: [no-new-privileges:true]`

### Production Path: Kubernetes

Docker Compose is suitable for single-node deployment (DGX Spark proof build). For Phase 2+ (departmental / multi-site), migration to Kubernetes provides:

- **Horizontal scaling:** Multiple replicas of stateless services (agent, embedding, FHIR publisher)
- **Service discovery:** DNS-based service resolution, load balancing
- **Rolling updates:** Zero-downtime deployment with readiness probes
- **Resource management:** CPU/memory requests and limits, GPU scheduling
- **Persistent storage:** StorageClass provisioners for NVMe, network storage
- **Secrets management:** Kubernetes Secrets for NGC_API_KEY, database credentials

### Discussion Questions — Chapter 5

1. Design a Kubernetes deployment for the agent that supports horizontal scaling of the inference pipeline while keeping the database as a StatefulSet.
2. What is the maximum theoretical container density on a DGX Spark given its 128 GB memory and 11 services?
3. How would you implement circuit-breaker patterns between the agent and NIM LLM service?

---

## Chapter 6 — PostgreSQL + pgvector

### Schema Design Principles

The agent's schema follows Third Normal Form (3NF) with selective denormalization:

- **studies** and **series** tables: normalized (series FK → studies)
- **findings** table: denormalized `details` JSONB column for workflow-specific data (avoids schema changes per workflow)
- **measurements** table: normalized (FK → findings), enabling cross-workflow measurement queries
- **embeddings** table: FK → studies, FK → findings (nullable), with vector(384) column

### Full DDL

```sql
CREATE EXTENSION IF NOT EXISTS vector;

CREATE TABLE studies (
    id                  SERIAL PRIMARY KEY,
    study_instance_uid  TEXT UNIQUE NOT NULL,
    patient_id          TEXT NOT NULL,
    patient_name        TEXT,
    study_date          DATE NOT NULL,
    study_description   TEXT,
    modality            TEXT NOT NULL,
    accession_number    TEXT,
    referring_physician TEXT,
    body_part           TEXT,
    num_series          INT DEFAULT 0,
    num_instances       INT DEFAULT 0,
    orthanc_id          TEXT,
    status              TEXT DEFAULT 'received',
    created_at          TIMESTAMPTZ DEFAULT NOW(),
    updated_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE series (
    id                  SERIAL PRIMARY KEY,
    series_instance_uid TEXT UNIQUE NOT NULL,
    study_id            INT REFERENCES studies(id) ON DELETE CASCADE,
    series_number       INT,
    series_description  TEXT,
    modality            TEXT NOT NULL,
    num_instances       INT DEFAULT 0,
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE findings (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    workflow         TEXT NOT NULL,
    finding_type    TEXT NOT NULL,
    finding_code    TEXT,
    location        TEXT,
    laterality      TEXT,
    severity        TEXT,
    confidence      FLOAT NOT NULL CHECK (confidence >= 0 AND confidence <= 1),
    is_positive     BOOLEAN DEFAULT TRUE,
    details         JSONB DEFAULT '{}',
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE measurements (
    id              SERIAL PRIMARY KEY,
    finding_id      INT REFERENCES findings(id) ON DELETE CASCADE,
    measurement_type TEXT NOT NULL,
    value           FLOAT NOT NULL,
    unit            TEXT NOT NULL,
    reference_range TEXT,
    flag            TEXT,
    prior_value     FLOAT,
    prior_date      DATE,
    delta_percent   FLOAT,
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE embeddings (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    finding_id      INT REFERENCES findings(id) ON DELETE SET NULL,
    level           TEXT NOT NULL,
    model_name      TEXT NOT NULL,
    embedding       vector(384) NOT NULL,
    metadata        JSONB DEFAULT '{}',
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE provenance (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    workflow         TEXT NOT NULL,
    model_id        TEXT NOT NULL,
    model_version   TEXT NOT NULL,
    model_arch      TEXT,
    inference_params JSONB DEFAULT '{}',
    input_uids      TEXT[] DEFAULT '{}',
    duration_ms     INT,
    gpu_memory_mb   INT,
    status          TEXT DEFAULT 'completed',
    error_message   TEXT,
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE worklist_entries (
    id              SERIAL PRIMARY KEY,
    study_id        INT REFERENCES studies(id) ON DELETE CASCADE,
    finding_id      INT REFERENCES findings(id) ON DELETE SET NULL,
    urgency         TEXT NOT NULL,
    priority        TEXT NOT NULL,
    notification    TEXT,
    routing         TEXT,
    acknowledged    BOOLEAN DEFAULT FALSE,
    acknowledged_by TEXT,
    acknowledged_at TIMESTAMPTZ,
    created_at      TIMESTAMPTZ DEFAULT NOW()
);
```

### pgvector Internals

pgvector stores vectors as fixed-length arrays of float4 (32-bit). The vector(384) column occupies 384 × 4 + 8 = 1,544 bytes per row (plus PostgreSQL tuple header overhead).

Distance operators:
- `<->` : L2 distance (Euclidean)
- `<=>` : Cosine distance (1 - cosine_similarity)
- `<#>` : Negative inner product

The agent uses cosine distance (`<=>`) for all similarity queries, consistent with BiomedCLIP's training objective.

### HNSW Index Tuning

```sql
CREATE INDEX idx_embeddings_hnsw ON embeddings
    USING hnsw (embedding vector_cosine_ops)
    WITH (m = 16, ef_construction = 64);
```

Parameters:
- **m = 16:** Maximum number of bi-directional links per node per layer. Higher m → better recall, larger index, slower inserts. Default is 16; range 4-64.
- **ef_construction = 64:** Size of the dynamic candidate list during index build. Higher → better index quality, slower build time. Must be ≥ 2m.
- **ef_search (runtime):** `SET hnsw.ef_search = 40;` Controls the size of the dynamic candidate list during search. Higher → better recall, slower queries. Default 40.

IVFFlat is the alternative: `CREATE INDEX ... USING ivfflat (embedding vector_cosine_ops) WITH (lists = 100);` and query with `SET ivfflat.probes = 10;`. Faster to build than HNSW but lower recall at the same query speed. HNSW is preferred for the agent's workload (moderate dataset size, high recall requirement).

### Hybrid Query Pattern

```sql
-- CTE filters by SQL predicates, then vector search on the filtered set
WITH candidates AS (
    SELECT e.study_id, e.embedding
    FROM embeddings e
    JOIN studies s ON e.study_id = s.id
    WHERE e.level = 'study'
      AND s.modality = 'CT'
      AND s.body_part = 'CHEST'
)
SELECT c.study_id, c.embedding <=> $1::vector AS distance
FROM candidates c
ORDER BY c.embedding <=> $1::vector
LIMIT 10;
```

PostgreSQL's query planner will use the HNSW index when the `ORDER BY ... LIMIT` pattern is present. For the CTE pattern, the planner may fall back to sequential scan on the filtered set — acceptable for small filtered sets but consider materialized views for large-scale deployments.

### Query Plan Analysis

```sql
EXPLAIN (ANALYZE, BUFFERS, FORMAT TEXT)
SELECT * FROM embeddings
ORDER BY embedding <=> '[0.1, 0.2, ...]'::vector
LIMIT 10;
```

Expected plan: `Index Scan using idx_embeddings_hnsw` with bounded heap scans. If the plan shows `Seq Scan`, verify the HNSW index exists and `hnsw.ef_search` is set.

### Discussion Questions — Chapter 6

1. What are the trade-offs of storing workflow-specific data in a JSONB column versus dedicated columns per workflow?
2. Calculate the approximate disk space required for 1 million embeddings at 384 dimensions.
3. Design a schema migration strategy that supports adding a new workflow without altering existing tables.

---

## Chapter 7 — Clinical Workflow Implementation

### MONAI Deploy MAP Lifecycle

A MAP (MONAI Deploy Application Package) follows this execution flow:

1. **Input:** DICOM files mounted at `/var/holoscan/input/`
2. **compose():** Defines the operator DAG (directed acyclic graph)
3. **Operators execute sequentially:** Each receives InputContext, produces OutputContext
4. **Output:** Results written to `/var/holoscan/output/`

```python
@resource(cpu=4, gpu=1, memory="16Gi")
class WorkflowApp(Application):
    def compose(self):
        self.add_flow(PreprocessOp(), InferenceOp())
        self.add_flow(InferenceOp(), PostprocessOp())
```

### CXR Rapid Findings — Implementation Detail

**Preprocessing:** Load CXR DICOM → resize to 224×224 (DenseNet-121 input) → normalize to [0,1] → ensure float32.

**DenseNet-121 architecture detail:** 4 dense blocks → global average pooling → FC(num_features, 5) with sigmoid activation (not softmax — multi-label, not multi-class). Each output is independent — a patient can have both pneumothorax AND pleural effusion.

**Threshold calibration:** Thresholds per class are set at the operating point that maximizes the Youden index (sensitivity + specificity - 1) on the validation set. Pneumothorax has a lower threshold (0.50) because the clinical cost of a false negative (missed collapsed lung) is higher than a false positive.

**GradCAM implementation:**
```python
cam = GradCAM(nn_module=model, target_layers="class_layers.relu")
heatmap = cam(x=input_tensor, class_idx=target_class)
```

GradCAM computes: (1) forward pass to get activations at target layer, (2) backward pass to get gradients of target class score w.r.t. activations, (3) global average pool the gradients, (4) weighted sum of activations, (5) ReLU to keep only positive contributions.

### CT Head Hemorrhage — Implementation Detail

**Two-stage pipeline:** Classification gate (DenseNet-121) → conditional segmentation (3D U-Net). The classification gate avoids running the expensive 3D segmentation on normal scans.

**Volume calculation:**
```python
volume_ml = np.sum(segmentation > 0) * np.prod(voxel_spacing) / 1000.0
```

Where `voxel_spacing` is in mm (from DICOM PixelSpacing and SpacingBetweenSlices), and division by 1000 converts mm³ to mL.

**Midline shift:** The agent computes the brain's anatomical midline (falx cerebri approximate position = axial center assuming RAS orientation), then measures the center of mass of the hemorrhage relative to this midline:

```python
shift_mm = abs((center_of_mass[0] - axial_center) * voxel_spacing[0])
```

**Brain Trauma Foundation decision tree:**
```python
if volume_ml > 30 or shift_mm > 5 or thickness_mm > 10:
    return ("critical", "P1")
elif volume_ml > 5:
    return ("urgent", "P2")
else:
    return ("routine", "P4")
```

### CT Chest Lung Nodule — Implementation Detail

**VDT derivation:** Assuming exponential growth: V(t) = V₀ × 2^(t/VDT). Solving for VDT:

```
V2/V1 = 2^(Δt/VDT)
ln(V2/V1) = (Δt/VDT) × ln(2)
VDT = (Δt × ln(2)) / ln(V2/V1)
```

VDT < 400 days is the commonly used threshold for suspicious growth rate (corresponds to approximately a 26% volume increase per year for a typical screening interval).

**Lung-RADS v2022 decision matrix:** The agent implements the full decision matrix including solid, ground-glass, and part-solid categories with different size thresholds for each. VDT < 400 upgrades the category.

**Genomics trigger:** Lung-RADS 4B or 4X → `check_genomics_trigger()` → POST to Nextflow API → Parabricks pipeline initiated.

### MRI Brain MS Lesion — Implementation Detail

**Preprocessing differences from CT:** MRI requires different normalization. Z-score normalization (`NormalizeIntensityd`) rather than HU windowing. N4 bias field correction (corrects spatial intensity inhomogeneity from RF coil sensitivity) is critical for accurate segmentation.

**Registration with ANTsPy:**
```python
result = ants.registration(fixed=current, moving=prior,
                           type_of_transform="SyNRA")
```

SyNRA = Rigid + Affine + SyN (symmetric normalization diffeomorphic registration). The diffeomorphic component handles non-linear brain deformation. Registration quality directly affects lesion matching accuracy.

**Lesion matching:** After warping prior lesion masks to current space, overlap (Dice coefficient) > 0.3 between a current and prior lesion indicates a match. New lesion = current lesion with no matching prior. Enlarging lesion = matched but volume increase > 20%.

### Discussion Questions — Chapter 7

1. Why is sigmoid activation used instead of softmax for CXR multi-label classification?
2. Derive the relationship between VDT and annual volume growth rate.
3. What registration failure modes could cause false-positive "new lesion" detections in MS tracking?

---

## Chapter 8 — LangGraph Agent Architecture

### StateGraph Internals

LangGraph's `StateGraph` compiles to a directed graph where:

- **State channels:** Each field in `AgentState` is a channel. Channels can have reducers that merge values from parallel branches.
- **Nodes:** Functions `f(state) → partial_state` that return updates to specific channels.
- **Edges:** Unconditional (always follow) or conditional (function evaluates state to choose next node).

```python
from langgraph.graph import StateGraph, END
from typing import TypedDict, Annotated

class AgentState(TypedDict):
    study_id: int
    findings: list[dict]
    prior_studies: list[dict]
    similar_cases: list[dict]
    severity: str
    report: str
    provenance: dict

graph = StateGraph(AgentState)
graph.add_node("triage", triage_node)
graph.add_node("longitudinal", longitudinal_node)
graph.add_node("population", population_node)
graph.add_node("report", report_node)

graph.set_entry_point("triage")
graph.add_conditional_edges("triage", route_by_severity,
    {"full_pipeline": "longitudinal", "brief_report": "report"})
graph.add_edge("longitudinal", "population")
graph.add_edge("population", "report")
graph.add_edge("report", END)

app = graph.compile()
```

### Conditional Routing Logic

```python
def route_by_severity(state: AgentState) -> str:
    if state["severity"] in ("critical", "urgent"):
        return "full_pipeline"
    return "brief_report"
```

### Tool Binding with MCP

Tools are declared with Pydantic schemas and bound via the Model Context Protocol:

```python
from langchain_core.tools import tool
from pydantic import BaseModel, Field

class SimilarStudyInput(BaseModel):
    embedding: list[float] = Field(description="384-dim query vector")
    modality: str = Field(description="CT, MR, CR, DX")
    limit: int = Field(default=10, description="Max results")

@tool(args_schema=SimilarStudyInput)
def search_similar_studies(embedding, modality, limit=10):
    """Search for studies with similar imaging characteristics."""
    return db.hybrid_similarity_search(embedding, modality, limit)
```

MCP publishes a tool manifest (name, description, parameters JSON schema) that any compatible agent runtime can discover and invoke.

### Checkpointing

LangGraph supports state checkpointing for:
- **Resumability:** If the agent fails mid-graph, resume from the last checkpoint
- **Human-in-the-loop:** Pause at a node, present state to user, await approval before continuing
- **Debugging:** Replay graph execution from any checkpoint

Backends: SQLite (proof build), PostgreSQL (production — reuses the agent's existing database).

### Agent Personas

Each agent node uses a tailored system prompt:

```python
TRIAGE_SYSTEM_PROMPT = """You are a triage radiologist assistant.
Your role: assess finding severity and determine routing.
Use Brain Trauma Foundation guidelines for hemorrhage.
Use Lung-RADS v2022 for lung nodules.
Output: severity classification and recommended analysis path."""
```

Different prompts for triage vs. longitudinal (compare to priors, calculate deltas) vs. population (interpret cohort data, outcome statistics) vs. report (synthesize all evidence into structured clinical narrative).

### Discussion Questions — Chapter 8

1. Design a fan-out/fan-in pattern where longitudinal and population nodes execute in parallel and their results are merged before the report node.
2. What are the trade-offs of using PostgreSQL vs. Redis as a LangGraph checkpoint backend?
3. How would you implement rate limiting for LLM tool calls to avoid overwhelming the NIM service?

---

## Chapter 9 — Embedding Models and Vector Retrieval

### Contrastive Learning and CLIP

CLIP (Contrastive Language-Image Pretraining) learns a joint embedding space for images and text. The training objective minimizes the InfoNCE loss:

```
L = -log(exp(sim(I_i, T_i)/τ) / Σ_j exp(sim(I_i, T_j)/τ))
```

Where sim() is cosine similarity between image embedding I and text embedding T, τ is a learnable temperature parameter, and the sum in the denominator runs over all image-text pairs in the batch. Matched pairs (I_i, T_i) are pulled together; mismatched pairs are pushed apart.

### BiomedCLIP Architecture

**Reference:** Zhang et al., "BiomedCLIP: A Multimodal Biomedical Foundation Model Pretrained from Fifteen Million Scientific Image-Text Pairs," 2023.

- **Text encoder:** PubMedBERT (BERT pretrained on PubMed abstracts)
- **Image encoder:** ViT-B/16 (Vision Transformer, base variant, 16×16 patch size)
- **Projection:** Both encoders output vectors projected to 384-dimensional shared space

The 384-dim output comes from the projection head design. The image encoder's CLS token (768-dim for ViT-B) is projected to 384-dim via a linear layer. This dimensionality balances representation capacity against storage and search efficiency.

### Vision Transformer (ViT) Detail

1. **Patch embedding:** Image (224×224) split into 14×14 = 196 patches of 16×16 pixels. Each patch flattened (16×16×channels = 768 for RGB) and linearly projected to embedding dimension d=768.
2. **Positional encoding:** Learnable 1D position embeddings added to patch embeddings.
3. **CLS token:** A learnable token prepended to the sequence. Its output serves as the image representation.
4. **Transformer blocks (12 for ViT-B):** Each block applies multi-head self-attention (MHSA) and feed-forward network (FFN) with layer normalization:

```
z'_l = MHSA(LN(z_{l-1})) + z_{l-1}
z_l = FFN(LN(z'_l)) + z'_l
```

**Scaled dot-product attention:**
```
Attention(Q, K, V) = softmax(QK^T / √d_k) × V
```

Where Q, K, V are query, key, value matrices derived from linear projections of the input. d_k is the key dimension. The scaling by √d_k prevents the dot products from growing large in magnitude, which would push the softmax into saturated regions.

### Three-Level Embedding Strategy

1. **Study-level:** Mean-pool all series-level embeddings for the study. Used for broad case-matching.
2. **Series-level:** For each series, sample representative slices, encode each, mean-pool. Used for protocol-specific matching.
3. **Lesion-level:** Crop ROI around individual findings, resize to 224×224, encode. Used for finding-specific matching ("nodules that look like this one").

### Vector Quantization

For scale (millions of embeddings), approximate nearest neighbor (ANN) techniques reduce storage and search cost:

**Product Quantization (PQ):** Split 384-dim vector into m subvectors (e.g., m=48, 8 dims each). Quantize each subvector to its nearest centroid from a learned codebook of k centroids (typically k=256). Store only the centroid indices (48 bytes vs. 1,536 bytes for the original). Distance approximated by summing precomputed sub-distances.

**HNSW + PQ:** Build the HNSW graph on quantized vectors for faster construction and smaller memory footprint, with optional reranking using exact vectors for the top candidates.

At the agent's expected scale (tens of thousands to low millions of embeddings), exact HNSW with full vectors is sufficient. PQ becomes relevant at Phase 3+ (multi-site, millions of studies).

### Embedding Drift

If the patient population or imaging protocols change over time, the embedding distribution shifts. Monitoring strategies:

- Track mean embedding per week/month — large shifts indicate distribution change
- Monitor retrieval recall on a held-out validation set
- Periodic reindexing when the distribution shift exceeds a threshold

### Discussion Questions — Chapter 9

1. Derive the storage and memory savings of PQ with m=48, k=256 compared to exact vectors at 384 dimensions.
2. Why does BiomedCLIP use 384-dim embeddings instead of the encoder's native 768 dimensions?
3. How would you implement an embedding update pipeline that recomputes embeddings when a new model version is deployed without downtime?

---

## Chapter 10 — Clinical Interoperability

### FHIR R4 Resource Model

The agent creates a FHIR Bundle containing:

```json
{
  "resourceType": "Bundle",
  "type": "transaction",
  "entry": [
    { "resource": { "resourceType": "DiagnosticReport", ... } },
    { "resource": { "resourceType": "Observation", ... } },
    { "resource": { "resourceType": "Observation", ... } }
  ]
}
```

**DiagnosticReport** structure:
- `status`: "final" (post-review) or "preliminary" (pre-review)
- `category`: LOINC 18726-0 (Radiology studies)
- `code`: LOINC code for the specific exam type
- `subject`: Reference to Patient resource
- `imagingStudy`: Reference to ImagingStudy (maps to DICOM Study UID)
- `result`: Array of references to Observation resources
- `conclusion`: Free-text summary (LLM-generated)
- `presentedForm`: Base64-encoded PDF/HTML report attachment

**Observation** structure:
- `code`: LOINC code for the measurement type
- `valueQuantity`: `{"value": 25.3, "unit": "mL", "system": "http://unitsofmeasure.org"}`
- `interpretation`: Coded clinical significance (e.g., high, critical)
- `component`: Multi-value observations (e.g., nodule with separate diameter, volume, VDT components)

### Terminology Binding

**SNOMED CT** hierarchy: Concepts have IS-A relationships (e.g., "Subdural hemorrhage" IS-A "Intracranial hemorrhage" IS-A "Hemorrhage"). The agent codes findings at the most specific level available.

**LOINC:** Structured along axes: Component/Property/Time/System/Scale/Method. Example: "Volume of hemorrhage in brain" = Component(hemorrhage volume) + System(brain) + Scale(quantitative).

### HL7v2 Integration (Legacy Pathway)

Many hospitals still run HL7v2 (pipe-delimited messages). Common messages:

- **ADT^A01 (Admit):** New patient registration → creates Patient resource
- **ORM^O01 (Order):** Imaging order placed → creates ServiceRequest
- **ORU^R01 (Result):** Radiology report → the agent's findings could be wrapped here

Integration engines (Mirth Connect, Rhapsody) translate between HL7v2 and FHIR. A FHIR façade pattern exposes a FHIR API while converting to/from HL7v2 behind the scenes.

### IHE AI Results Profile

The Integrating the Healthcare Enterprise (IHE) organization publishes integration profiles. The IHE-AI Results supplement (under development) defines how AI findings should be communicated to PACS, EHR, and worklist systems — including provenance, confidence scores, and links to evidence images. The agent's architecture aligns with this profile.

### SMART on FHIR Security

SMART (Substitutable Medical Applications, Reusable Technologies) on FHIR provides:
- OAuth2 authorization with FHIR-specific scopes (e.g., `patient/Observation.read`)
- Launch context (which patient, which encounter)
- Token-based access control

The agent's FHIR publisher authenticates with the EHR's FHIR server using client credentials (machine-to-machine OAuth2 flow).

### Discussion Questions — Chapter 10

1. Design the FHIR Observation resources needed to represent a lung nodule with diameter, volume, VDT, and Lung-RADS score as components.
2. How would you implement a FHIR façade over the agent's PostgreSQL database to enable external EHR systems to query findings?
3. What IHE profiles are relevant for automated AI-to-PACS worklist prioritization?

---

## Chapter 11 — RAG Architecture and LLM Serving

### Transformer Architecture

The transformer (Vaswani et al., 2017) processes sequences using self-attention:

```
MultiHead(Q,K,V) = Concat(head_1, ..., head_h) × W^O
where head_i = Attention(QW_i^Q, KW_i^K, VW_i^V)
```

Each transformer block: LayerNorm → MHSA → Residual → LayerNorm → FFN → Residual.

### Llama 3 8B Architecture

- 32 transformer layers
- 32 attention heads, 8 KV heads (Grouped Query Attention — GQA)
- Hidden dimension: 4096
- FFN dimension: 14336 (SwiGLU activation)
- RoPE (Rotary Position Embedding) for positional encoding
- Context window: 8192 tokens
- Vocabulary: 128K tokens (byte-level BPE)

**GQA:** Instead of separate KV projections per head, GQA groups 4 query heads per KV head (32 Q / 8 KV). This reduces KV cache memory by 4× with minimal quality loss — critical for inference efficiency.

**RoPE:** Encodes position through rotation matrices applied to Q and K vectors. The rotation angle is proportional to position, enabling relative position encoding that generalizes to unseen sequence lengths.

### NVIDIA NIM Serving

NIM handles model serving with:
- **Model loading:** Loads weights into GPU memory, initializes KV cache
- **Continuous batching:** Dynamically adds/removes requests from the batch as they arrive/complete — maximizes GPU utilization vs. static batching
- **KV cache management:** Paged attention (inspired by vLLM) — KV cache allocated in pages, enabling non-contiguous memory allocation and reducing fragmentation
- **Speculative decoding:** Small draft model generates candidate tokens; large model verifies in parallel — can accelerate generation 2-3× when draft acceptance rate is high

NIM container on DGX Spark:
```yaml
image: nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark
```

The `-dgx-spark` tag indicates ARM64 optimization for Grace-Blackwell.

### RAG Pipeline Design

The agent's RAG pipeline follows a retriever → generator pattern (no explicit reranker for latency reasons):

1. **Dense retrieval:** Query embeddings search pgvector for similar studies
2. **SQL retrieval:** Direct database queries for structured findings, measurements, provenance
3. **Evidence assembly:** All retrieved evidence formatted into a structured prompt
4. **Generation:** NIM-served Llama 3 generates the report with evidence citations

**Context window budget (8192 tokens):**
- System prompt: ~500 tokens
- Current findings: ~300-500 tokens
- Prior measurements: ~200-400 tokens
- Similar cases: ~500-1000 tokens (top 5 cases, abbreviated)
- Guidelines: ~500-800 tokens
- Generation headroom: ~5000 tokens for output

### Prompt Engineering

System message template:
```
You are a radiology AI assistant generating clinical reports.
RULES:
1. Only include facts supported by the evidence below.
2. Never add information not present in the evidence.
3. Cite specific measurements with units.
4. Flag critical findings prominently.
5. Compare to prior studies when available.
6. Reference applicable clinical guidelines.

EVIDENCE:
{evidence_block}

Generate a structured clinical report with sections:
1. Clinical Indication
2. Findings
3. Measurements
4. Comparison to Prior
5. Impression
6. Recommendation
```

### Evaluation with RAGAS

RAGAS (Retrieval Augmented Generation Assessment) provides metrics:
- **Faithfulness:** What fraction of claims in the generated report are supported by the retrieved context? (Ground truth: evidence passages)
- **Answer relevancy:** Does the report address the clinical question?
- **Context recall:** Did the retriever find all relevant evidence?

### Discussion Questions — Chapter 11

1. Calculate the KV cache memory requirement for Llama 3 8B at full context length (8192 tokens) with GQA.
2. Design a hybrid retrieval strategy combining BM25 keyword search with dense embedding search using Reciprocal Rank Fusion.
3. How would you detect and prevent hallucinated measurements in the generated report?

---

## Chapter 12 — Observability Stack

### Prometheus Architecture

Prometheus operates a **pull model** — it scrapes HTTP endpoints at configured intervals:

```yaml
# prometheus.yml
scrape_configs:
  - job_name: 'imaging-agent'
    static_configs:
      - targets: ['agent:8000']
    scrape_interval: 15s
  - job_name: 'dcgm'
    static_configs:
      - targets: ['dcgm-exporter:9400']
    scrape_interval: 15s
```

Time series database (TSDB) stores samples as (timestamp, value) pairs, indexed by metric name and label set. Retention configured at 30 days (`--storage.tsdb.retention.time=30d`).

### Metric Types

| Type | Description | Example |
|---|---|---|
| Counter | Monotonically increasing | `studies_processed_total` |
| Gauge | Arbitrary value (up/down) | `gpu_utilization_percent` |
| Histogram | Distribution in configurable buckets | `inference_duration_seconds` |
| Summary | Distribution with configurable quantiles | `request_duration_seconds` |

Instrumentation in Python:
```python
from prometheus_client import Counter, Histogram, Gauge

STUDIES_TOTAL = Counter('studies_processed_total', 'Total studies processed', ['workflow'])
INFERENCE_DURATION = Histogram('inference_duration_seconds', 'Inference time',
                               ['workflow'], buckets=[0.5, 1, 2, 5, 10, 30, 60, 120])
GPU_UTIL = Gauge('gpu_utilization_percent', 'Current GPU utilization')
```

### PromQL Queries

Key queries for the agent:

```promql
# Per-workflow inference latency (p95 over 5 minutes)
histogram_quantile(0.95, rate(inference_duration_seconds_bucket[5m]))

# Studies processed per hour
rate(studies_processed_total[1h]) * 3600

# GPU utilization (from DCGM)
DCGM_FI_DEV_GPU_UTIL{gpu="0"}
```

### DCGM Telemetry Fields

| Field | Description |
|---|---|
| DCGM_FI_DEV_GPU_UTIL | GPU compute utilization (%) |
| DCGM_FI_DEV_MEM_COPY_UTIL | Memory controller utilization (%) |
| DCGM_FI_DEV_GPU_TEMP | GPU temperature (°C) |
| DCGM_FI_DEV_POWER_USAGE | Power consumption (W) |
| DCGM_FI_DEV_FB_FREE | Free framebuffer memory (MiB) |
| DCGM_FI_DEV_FB_USED | Used framebuffer memory (MiB) |
| DCGM_FI_DEV_ECC_SBE_VOL | Single-bit ECC errors (volatile) |
| DCGM_FI_DEV_XID_ERRORS | XID error count |

### OpenTelemetry Integration

For distributed tracing across containers, OpenTelemetry provides:
- **Spans:** Timed operations (e.g., "dicom_listener.process_study", "agent.triage_node")
- **Trace context propagation:** W3C traceparent header passed between HTTP calls
- **Exporters:** Send traces to Jaeger, Zipkin, or Grafana Tempo

This enables end-to-end latency breakdown: how long did DICOM ingestion take? Inference? Database queries? Report generation?

### SLA Definition

| Workflow | p50 Target | p95 Target | p99 Target |
|---|---|---|---|
| CXR Rapid Findings | < 10s | < 20s | < 30s |
| CT Head Hemorrhage | < 45s | < 75s | < 90s |
| CT Chest Lung Nodule | < 2min | < 4min | < 5min |
| MRI Brain MS Lesion | < 2min | < 4min | < 5min |

### Discussion Questions — Chapter 12

1. Design a Grafana alerting rule that fires when the p95 inference latency exceeds the SLA target for any workflow.
2. What DCGM metrics would you monitor to predict GPU hardware failure before it occurs?
3. How would you implement distributed tracing across the agent's 11 containers without modifying application code?

---

## Chapter 13 — HCLS AI Factory

### Multi-Agent Coordination

The AI Factory agents communicate through:
- **Shared PostgreSQL database:** Common schema for cross-agent queries
- **FHIR ServiceRequest:** Trigger messages between agents (imaging → genomics)
- **Event bus:** Study completion events, finding alerts, genomics results (Kafka or Redis Streams at scale)

### Parabricks Pipeline

NVIDIA Parabricks GPU-accelerates genomics:
- **BWA-MEM2:** Read alignment to reference genome — GPU kernels for Smith-Waterman alignment
- **HaplotypeCaller / DeepVariant:** Variant calling (identify SNPs, indels, structural variants)
- **Performance:** 30× whole-genome sequencing in ~10 minutes on 8×A100 GPUs (vs. ~30 hours CPU)

The imaging-to-genomics trigger:
```
Lung nodule Lung-RADS 4B+ → FHIR ServiceRequest → Nextflow pipeline → Parabricks → VCF output → Variant annotation
```

### Federated Learning: FedAvg Algorithm

```
For each communication round t:
    Server sends global model w_t to selected clients
    Each client k:
        w_k^{t+1} = w_t - η∇L_k(w_t)  (local training for E epochs)
    Server aggregates:
        w_{t+1} = Σ_k (n_k/n) × w_k^{t+1}  (weighted average by local dataset size)
```

**NVIDIA FLARE** implements this with:
- **Scatter-and-Gather workflow:** Server scatters model to clients, clients train locally, server gathers and aggregates
- **Privacy modules:** Differential privacy (add calibrated noise to gradients), secure aggregation (cryptographic protocols preventing server from seeing individual updates)
- **Model filters:** Client-side filters that clip gradients, quantize updates, or apply sparsification before transmission

### RECIST 1.1 Criteria

Response Evaluation Criteria in Solid Tumors:
- **Complete Response (CR):** Disappearance of all target lesions
- **Partial Response (PR):** ≥30% decrease in sum of longest diameters
- **Progressive Disease (PD):** ≥20% increase in sum of longest diameters AND ≥5mm absolute increase
- **Stable Disease (SD):** Neither PR nor PD

The imaging agent's automated measurements can serve as RECIST endpoints in clinical trials — replacing manual radiologist measurements with reproducible, automated quantification.

### Nextflow DSL2 Pipeline Orchestration

```groovy
// main.nf
nextflow.enable.dsl = 2

include { CT_HEAD_HEMORRHAGE } from './modules/ct_head_hemorrhage'
include { CT_CHEST_LUNG_NODULE } from './modules/ct_chest_lung_nodule'
include { CXR_RAPID_FINDINGS } from './modules/cxr_rapid_findings'
include { MRI_BRAIN_MS_LESION } from './modules/mri_brain_ms_lesion'

workflow {
    study_ch = Channel.fromPath(params.input_dir)

    // Route by modality and body part
    study_ch.branch {
        ct_head: it.modality == 'CT' && it.body_part == 'HEAD'
        ct_chest: it.modality == 'CT' && it.body_part == 'CHEST'
        cxr: it.modality in ['CR', 'DX']
        mri_brain: it.modality == 'MR' && it.body_part == 'BRAIN'
    }.set { routed }

    CT_HEAD_HEMORRHAGE(routed.ct_head)
    CT_CHEST_LUNG_NODULE(routed.ct_chest)
    CXR_RAPID_FINDINGS(routed.cxr)
    MRI_BRAIN_MS_LESION(routed.mri_brain)
}
```

### Scaling Analysis

| GPU Tier | Studies/Hour (CXR) | Studies/Hour (CT Head) | Concurrent Pipelines |
|---|---|---|---|
| DGX Spark (1 GPU) | ~120 | ~40 | 1 |
| DGX B200 (8 GPU) | ~960 | ~320 | 8 |
| DGX SuperPOD (256 GPU) | ~30,000 | ~10,000 | 256 |

Estimates assume single-study-per-GPU sequential processing. Throughput increases linearly with GPU count for embarrassingly parallel study processing.

### Discussion Questions — Chapter 13

1. Design a differential privacy mechanism for the hemorrhage detection model using NVIDIA FLARE that provides ε=1 privacy guarantee.
2. How would you implement RECIST 1.1 automated measurement with the lung nodule workflow?
3. Calculate the total compute cost (hardware + NVAIE license) for a Phase 3 multi-site deployment processing 500 studies/hour across 4 hospitals.

---

## Chapter 14 — Regulatory, Safety, and Deployment Engineering

### FDA SaMD Risk Categorization

The FDA categorizes SaMD based on two axes:

| | Treat/Diagnose | Drive Management | Inform |
|---|---|---|---|
| Critical | IV (highest) | III | II |
| Serious | III | II | I (lowest) |
| Non-serious | II | I | I |

The imaging agent as triage decision support for serious conditions (hemorrhage, cancer) falls in Category II-III depending on the specific claim.

### Predetermined Change Control Plan (PCCP)

A PCCP consists of:

**SaMD Pre-Specifications (SPS):** Describe the types of changes anticipated:
- "The model may be retrained on expanded datasets including new scanner manufacturers"
- "Confidence thresholds may be adjusted based on real-world performance data"
- "New finding types may be added to the CXR workflow"

**Algorithm Change Protocol (ACP):** Define the validation methodology:
- "Any retrained model must achieve sensitivity ≥95% and specificity ≥90% on the fixed validation set"
- "Threshold changes must be validated on a held-out calibration set with ≥1000 positive cases"
- "Changes exceeding the SPS scope require a new regulatory submission"

### Clinical Validation Methodology

| Metric | Formula | Interpretation |
|---|---|---|
| Sensitivity | TP / (TP + FN) | Fraction of true positives detected |
| Specificity | TN / (TN + FP) | Fraction of true negatives correctly identified |
| PPV | TP / (TP + FP) | Fraction of positive predictions that are correct |
| NPV | TN / (TN + FN) | Fraction of negative predictions that are correct |
| AUC-ROC | Area under ROC curve | Overall discrimination ability |
| AUC-PR | Area under precision-recall curve | Performance on imbalanced datasets (preferred for rare findings) |

AUC-PR is preferred over AUC-ROC for medical imaging because class imbalance is extreme (e.g., <1% of CXRs have pneumothorax). AUC-ROC can appear favorable even with many false positives when the negative class is large.

### Bias Assessment

Demographic subgroup analysis: compute all metrics stratified by age group, sex, race/ethnicity, scanner manufacturer, and site. Report performance gaps with 95% confidence intervals.

Fairness metrics:
- **Equalized odds:** P(Ŷ=1|Y=1, A=a) = P(Ŷ=1|Y=1, A=b) for all groups a, b (equal sensitivity)
- **Demographic parity:** P(Ŷ=1|A=a) = P(Ŷ=1|A=b) (equal positive prediction rate)
- **Calibration:** P(Y=1|Ŷ=p) = p for all groups (confidence scores are accurate across groups)

### IEC 62304 Software Lifecycle

Medical device software development standard requiring:
- **Software development plan** (architecture, design, testing strategy)
- **Software requirements specification** (functional, performance, safety)
- **Software architecture design** (modules, interfaces, data flow)
- **Unit testing, integration testing, system testing** (with traceability to requirements)
- **Software maintenance plan** (bug fixes, updates, monitoring)

Classification: Class A (no injury), B (non-serious injury), C (death or serious injury). The agent's triage function is Class C.

### Deployment Patterns

**Shadow mode (parallel run):** New model processes same studies as production model. Results compared offline. No impact on clinical workflow. Validates performance before promotion.

**Canary deployment:** New model serves a small percentage of studies (e.g., 5%). Monitor error rates and latency. Gradually increase if metrics are stable.

**Blue-green deployment:** Two identical environments. Switch traffic from blue (current) to green (new) atomically. Instant rollback by switching back.

### Responsible AI

**Explainability methods:**
- GradCAM: Where is the model looking? (spatial explanation)
- SHAP: Which input features contribute most to the prediction? (feature attribution)
- Integrated Gradients: Axiomatic attribution method — satisfies sensitivity and implementation invariance

**Uncertainty quantification:**
- MC Dropout: Run inference N times with different dropout masks, compute prediction variance
- Deep Ensemble: Train N independent models, compute prediction variance across ensemble
- Epistemic uncertainty (model uncertainty) vs. aleatoric uncertainty (data noise)

### Discussion Questions — Chapter 14

1. Draft a PCCP SPS and ACP for the CT head hemorrhage workflow that allows threshold adjustment and model retraining without a new 510(k).
2. How would you design a demographic bias audit for the CXR workflow given that training data underrepresents certain age groups?
3. Compare MC Dropout vs. Deep Ensemble for uncertainty quantification in terms of inference cost on the DGX Spark.

---

## Discussion Questions — Comprehensive

These 20 questions require synthesis across multiple chapters.

1. Design an end-to-end data flow from DICOM C-STORE receipt through LangGraph agent processing to FHIR DiagnosticReport output. Identify every container involved and every database write.

2. Calculate the total GPU memory budget for running the hemorrhage 3D U-Net, the NIM LLM (Llama 3 8B), and the BiomedCLIP embedding model concurrently on the DGX Spark's 128 GB unified memory.

3. Design a federated learning protocol using NVIDIA FLARE that improves hemorrhage detection across 5 hospitals while preserving ε=1 differential privacy. What FLARE components would you configure?

4. How would you modify the Lung-RADS classification algorithm to handle the case where prior study registration fails? What fallback strategy preserves clinical safety?

5. Propose a schema migration that adds a new "cardiac_function" workflow without modifying existing tables or breaking existing queries.

6. Design a Grafana dashboard for the hemorrhage workflow showing: p95 latency, throughput, CRITICAL finding rate, GPU utilization, and NIM token throughput.

7. Compare the annotation cost of training data for classification (DenseNet-121) vs. segmentation (3D U-Net) vs. detection (RetinaNet). How does this affect model iteration speed?

8. Design an A/B test protocol that compares radiologist performance with and without the agent's triage assistance while maintaining patient safety.

9. How would you implement a "second opinion" pattern where two independent models must agree on a CRITICAL classification before triggering an alert?

10. Calculate the cosine similarity search latency for 1 million embeddings using HNSW with m=16, ef_search=100. Compare to IVFFlat with lists=1000, probes=20.

11. Design the FHIR resources and terminology bindings needed to represent a longitudinal MS lesion tracking report with new, enlarging, and stable lesion counts.

12. How would you extend the RAG pipeline to include genomic variant data from Parabricks in the evidence context for a lung cancer staging report?

13. Propose a TensorRT optimization strategy for the agent's models. Which models benefit most from INT8 quantization? Which require FP16 minimum?

14. Design a container health check cascade that detects and automatically recovers from a NIM LLM out-of-memory error without losing in-flight studies.

15. How would the agent's architecture change if deployed on a DGX B200 (8 GPUs) instead of DGX Spark (1 GPU)? What parallelism strategies become available?

16. Design a provenance query that reconstructs the complete processing history for a specific patient across all modalities and workflows over 2 years.

17. How would you implement real-time embedding drift detection and automated reindexing without disrupting similarity search availability?

18. Compare the security implications of SMART on FHIR vs. API key authentication for the agent's FHIR publisher in a multi-tenant hospital deployment.

19. Design a Nextflow pipeline that implements the cross-modal trigger: imaging Lung-RADS 4B+ → genomics → combined report. Include error handling for Parabricks pipeline failure.

20. Estimate the total cost of ownership (hardware + NVAIE + operational) for a Phase 2 departmental deployment processing 200 studies/hour across all 4 workflows.

---

## References

### Deep Learning Architectures
- Huang et al., "Densely Connected Convolutional Networks," CVPR 2017
- Ronneberger et al., "U-Net: Convolutional Networks for Biomedical Image Segmentation," MICCAI 2015
- Çiçek et al., "3D U-Net: Learning Dense Volumetric Segmentation from Sparse Annotation," MICCAI 2016
- Myronenko, "3D MRI Brain Tumor Segmentation Using Autoencoder Regularization," BrainLes@MICCAI 2018
- Lin et al., "Focal Loss for Dense Object Detection," ICCV 2017
- Lin et al., "Feature Pyramid Networks for Object Detection," CVPR 2017
- Hatamizadeh et al., "Swin UNETR: Swin Transformers for Semantic Segmentation of Brain Tumors," BrainLes@MICCAI 2021

### Foundation Models
- Radford et al., "Learning Transferable Visual Models From Natural Language Supervision," ICML 2021 (CLIP)
- Zhang et al., "BiomedCLIP: A Multimodal Biomedical Foundation Model," 2023
- Dosovitskiy et al., "An Image is Worth 16x16 Words: Transformers for Image Recognition at Scale," ICLR 2021 (ViT)
- Vaswani et al., "Attention Is All You Need," NeurIPS 2017
- Touvron et al., "Llama 2: Open Foundation and Fine-Tuned Chat Models," 2023
- Meta AI, "Llama 3 Model Card," 2024

### Medical Imaging Standards
- DICOM Standard PS3.3: Information Object Definitions
- DICOM Standard PS3.4: Service Class Specifications
- DICOM Standard PS3.5: Data Structures and Encoding
- DICOM Standard PS3.18: Web Services (DICOMweb)
- HL7 FHIR R4 Specification (hl7.org/fhir/R4)
- IHE Radiology Technical Framework

### Clinical Guidelines
- ACR Lung-RADS v2022 Assessment Categories
- Brain Trauma Foundation, "Guidelines for the Management of Severe Traumatic Brain Injury," 4th Edition
- Thompson et al., "Diagnosis of Multiple Sclerosis: 2017 Revisions of the McDonald Criteria," Lancet Neurology 2018
- RECIST 1.1: Eisenhauer et al., "New Response Evaluation Criteria in Solid Tumours," European Journal of Cancer 2009

### Regulatory
- FDA, "Artificial Intelligence/Machine Learning (AI/ML)-Based Software as a Medical Device (SaMD) Action Plan," 2021
- FDA, "Marketing Submission Recommendations for a Predetermined Change Control Plan for AI/ML-Enabled Device Software Functions," 2023
- IEC 62304: Medical Device Software — Software Life Cycle Processes
- ISO 14971: Medical Devices — Application of Risk Management
- IEC 62366: Medical Devices — Usability Engineering

### NVIDIA Platforms
- NVIDIA DGX Spark Technical Specifications
- NVIDIA MONAI Deploy SDK Documentation
- NVIDIA NIM Documentation
- NVIDIA Parabricks Documentation
- NVIDIA FLARE Documentation (Apache 2.0)
- NVIDIA BioNeMo Service Documentation

---

> *This document was created for the HCLS AI Factory — Imaging Intelligence Agent.*
> *Apache 2.0 License | Author: Adam Jones | February 2026*
