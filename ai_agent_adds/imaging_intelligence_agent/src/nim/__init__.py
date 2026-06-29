"""NIM client layer for the Imaging Intelligence Agent.

Provides clients for four NVIDIA NIM microservices:
  - VISTA-3D: 3D medical image segmentation (127+ anatomical classes)
  - MAISI: Synthetic CT volume generation
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
from .service_manager import NIMServiceManager
from .vilam3_client import VILAM3Client
from .vista3d_client import VISTA3DClient

__all__ = [
    "BaseNIMClient",
    "VISTA3DClient",
    "MAISIClient",
    "VILAM3Client",
    "LlamaLLMClient",
    "NIMServiceManager",
]
