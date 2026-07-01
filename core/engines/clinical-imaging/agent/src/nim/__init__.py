"""NIM client layer for the Imaging Intelligence Agent.

Provides clients for seven NVIDIA NIM microservices:
  - VISTA-3D: 3D medical image segmentation (127+ anatomical classes)
  - NV-Segment-CT: 3D CT segmentation with 132 classes incl. 7 tumors
  - MAISI: Synthetic CT volume generation
  - NV-Generate-CT: Synthetic CT volume generation (rectified flow, 33x faster)
  - NV-Generate-MR: Synthetic MRI volume generation (multi-contrast)
  - VILA-M3: Visual language model for medical image understanding
  - LLM: Llama-3 text generation with Anthropic Claude fallback

All clients inherit from BaseNIMClient which provides:
  - Cached health checks with configurable interval
  - Exponential-backoff retry via tenacity
  - Automatic mock fallback when services are unavailable
"""

from .base import BaseNIMClient
from .llm_client import LlamaLLMClient
from .maisi_client import MAISIClient
from .nemotron_nano_client import NemotronNanoClient
from .nv_generate_ct_client import NVGenerateCTClient
from .nv_generate_mr_client import NVGenerateMRClient
from .nv_segment_ct_client import NVSegmentCTClient
from .service_manager import NIMServiceManager
from .vilam3_client import VILAM3Client
from .vista3d_client import VISTA3DClient

__all__ = [
    "BaseNIMClient",
    "VISTA3DClient",
    "NVSegmentCTClient",
    "MAISIClient",
    "NVGenerateCTClient",
    "NVGenerateMRClient",
    "VILAM3Client",
    "LlamaLLMClient",
    "NemotronNanoClient",
    "NIMServiceManager",
]
