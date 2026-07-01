# Imaging Intelligence Agent — Technology Licenses

Deployment tier license inventory for every model, framework, and service
used by the Imaging Intelligence Agent.

## Community Edition (20 technologies — all commercially usable)

| # | Technology | License | Source | Commercial Use |
|---|-----------|---------|--------|----------------|
| 1 | MONAI Core | Apache 2.0 | [GitHub](https://github.com/Project-MONAI/MONAI) | Yes |
| 2 | MONAI Label | Apache 2.0 | [GitHub](https://github.com/Project-MONAI/MONAILabel) | Yes |
| 3 | MONAI Deploy | Apache 2.0 | [GitHub](https://github.com/Project-MONAI/monai-deploy) | Yes |
| 4 | NV-Segment-CT | Open Model | [NGC](https://catalog.ngc.nvidia.com/orgs/nim/teams/nvidia/models/nv-segment-ct) | Yes |
| 5 | TorchXRayVision (DenseNet-121) | Apache 2.0 | [GitHub](https://github.com/mlmed/torchxrayvision) | Yes |
| 6 | FLARE (federated learning) | Apache 2.0 | [GitHub](https://github.com/NVIDIA/NVFlare) | Yes |
| 7 | Holoscan SDK | Apache 2.0 | [GitHub](https://github.com/nvidia-holoscan/holoscan-sdk) | Yes |
| 8 | BGE-small-en-v1.5 | MIT | [HuggingFace](https://huggingface.co/BAAI/bge-small-en-v1.5) | Yes |
| 9 | Milvus | Apache 2.0 | [GitHub](https://github.com/milvus-io/milvus) | Yes |
| 10 | Orthanc DICOM | GPLv3 | [Orthanc](https://www.orthanc-server.com/) | Yes (server) |
| 11 | OHIF Viewer | MIT | [GitHub](https://github.com/OHIF/Viewers) | Yes |
| 12 | PyDICOM | MIT | [PyPI](https://pypi.org/project/pydicom/) | Yes |
| 13 | highdicom | MIT | [PyPI](https://pypi.org/project/highdicom/) | Yes |
| 14 | SimpleITK | Apache 2.0 | [PyPI](https://pypi.org/project/SimpleITK/) | Yes |
| 15 | NiBabel | MIT | [PyPI](https://pypi.org/project/nibabel/) | Yes |
| 16 | PyRadiomics | BSD 3-Clause | [GitHub](https://github.com/AIM-Harvard/pyradiomics) | Yes |
| 17 | RDKit | BSD 3-Clause | [GitHub](https://github.com/rdkit/rdkit) | Yes |
| 18 | Streamlit | Apache 2.0 | [PyPI](https://pypi.org/project/streamlit/) | Yes |
| 19 | FastAPI | MIT | [PyPI](https://pypi.org/project/fastapi/) | Yes |
| 20 | Meta Llama 3 8B Instruct | Llama 3 Community | [NGC](https://catalog.ngc.nvidia.com/orgs/nvidia/teams/nim/containers/meta-llama3-8b-instruct) | Yes |

## Enterprise Edition (adds 8 technologies — requires NVIDIA AI Enterprise)

| # | Technology | License | Source | Commercial Use |
|---|-----------|---------|--------|----------------|
| 21 | NVIDIA VISTA-3D NIM | AI Enterprise | [NGC](https://catalog.ngc.nvidia.com/orgs/nvidia/teams/nim/containers/vista3d) | Yes (with license) |
| 22 | NVIDIA MAISI NIM | AI Enterprise | [NGC](https://catalog.ngc.nvidia.com/orgs/nvidia/teams/nim/containers/maisi) | Yes (with license) |
| 23 | VILA-M3 NIM | AI Enterprise | [NGC](https://catalog.ngc.nvidia.com/orgs/nvidia/teams/nim/containers/vilam3) | Yes (with license) |
| 24 | NeMo Retriever Embedding NIM | AI Enterprise | [NGC](https://catalog.ngc.nvidia.com/orgs/nvidia/teams/nim/models/nemo-retriever) | Yes (with license) |
| 25 | NeMo Guardrails | Apache 2.0 | [GitHub](https://github.com/NVIDIA/NeMo-Guardrails) | Yes |
| 26 | the AI platform | Commercial | [the platform](https://vastdata.com/platform) | Yes (with license) |
| 27 | NVIDIA Triton Inference Server | BSD 3-Clause | [GitHub](https://github.com/triton-inference-server/server) | Yes |
| 28 | NVIDIA DCGM | Apache 2.0 | [GitHub](https://github.com/NVIDIA/DCGM) | Yes |

## Research Edition (adds 2 technologies — noncommercial only)

| # | Technology | License | Source | Commercial Use |
|---|-----------|---------|--------|----------------|
| 29 | NV-Reason-CXR-3B | Noncommercial | [NGC](https://catalog.ngc.nvidia.com/orgs/nim/teams/nvidia/models/nv-reason-cxr-3b) | No |
| 30 | NV-Segment-CTMR | Noncommercial | [NGC](https://catalog.ngc.nvidia.com/orgs/nim/teams/nvidia/models/nv-segment-ctmr) | No |

## Tier Summary

| Tier | Total Technologies | GPU Required | Commercial Use | Key Additions |
|------|--------------------|-------------|----------------|---------------|
| Community | 20 | Yes (1x) | Yes | All open-source / open-model |
| Enterprise | 28 | Yes (1x+) | Yes (licensed) | VISTA-3D, MAISI, NeMo, the AI platform |
| Research | 30 | Yes (1x+) | No (restricted) | NV-Reason-CXR-3B, NV-Segment-CTMR |

## License References

- Apache 2.0: <https://www.apache.org/licenses/LICENSE-2.0>
- MIT: <https://opensource.org/licenses/MIT>
- BSD 3-Clause: <https://opensource.org/licenses/BSD-3-Clause>
- GPLv3: <https://www.gnu.org/licenses/gpl-3.0.html>
- Llama 3 Community License: <https://llama.meta.com/llama3/license/>
- NVIDIA AI Enterprise EULA: <https://www.nvidia.com/en-us/data-center/products/ai-enterprise/eula/>
- NVIDIA Open Model License: <https://developer.nvidia.com/open-model-license>
- NVIDIA Noncommercial License: <https://developer.nvidia.com/noncommercial-license>
