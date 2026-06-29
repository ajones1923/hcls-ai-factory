"""NIM Service Manager for the Imaging Intelligence Agent.

Centralizes lifecycle management, health monitoring, and access
to all four NIM service clients (VISTA-3D, MAISI, VILA-M3, LLM).

Supports three modes per service:
  - "available": local NIM endpoint is live
  - "cloud": NVIDIA cloud NIM endpoint is live (LLM and VLM only)
  - "mock": using mock/demo responses
  - "unavailable": no backend and mock disabled
"""

from typing import Dict, List, Optional

from loguru import logger

from .base import BaseNIMClient
from .llm_client import LlamaLLMClient
from .maisi_client import MAISIClient
from .vilam3_client import VILAM3Client
from .vista3d_client import VISTA3DClient


class NIMServiceManager:
    """Manages all NIM service clients for the Imaging Intelligence Agent.

    Creates and configures clients from ImagingSettings, provides
    aggregated health checks, and exposes typed properties for
    convenient access.

    Usage:
        from config.settings import settings
        manager = NIMServiceManager(settings)
        status = manager.check_all_services()
        result = manager.vista3d.segment("volume.nii.gz")
    """

    def __init__(self, settings, nvidia_api_key: Optional[str] = None):
        """Initialize all NIM clients from settings.

        Args:
            settings: ImagingSettings instance with NIM_*_URL and
                NIM_ALLOW_MOCK_FALLBACK configuration.
            nvidia_api_key: Optional NVIDIA API key for cloud NIM endpoints.
                If not provided, reads from settings.NVIDIA_API_KEY.
        """
        mock_enabled = getattr(settings, "NIM_ALLOW_MOCK_FALLBACK", True)

        # Resolve NVIDIA API key: explicit param > settings > None
        resolved_nvidia_key = nvidia_api_key or getattr(settings, "NVIDIA_API_KEY", None)
        cloud_url = getattr(settings, "NIM_CLOUD_URL", "https://integrate.api.nvidia.com/v1")
        cloud_llm_model = getattr(settings, "NIM_CLOUD_LLM_MODEL", "meta/llama-3.1-8b-instruct")
        cloud_vlm_model = getattr(settings, "NIM_CLOUD_VLM_MODEL", "meta/llama-3.2-11b-vision-instruct")

        self._vista3d = VISTA3DClient(
            base_url=settings.NIM_VISTA3D_URL,
            mock_enabled=mock_enabled,
        )
        self._maisi = MAISIClient(
            base_url=settings.NIM_MAISI_URL,
            mock_enabled=mock_enabled,
        )
        self._vilam3 = VILAM3Client(
            base_url=settings.NIM_VILAM3_URL,
            mock_enabled=mock_enabled,
            nvidia_api_key=resolved_nvidia_key,
            cloud_url=cloud_url,
            cloud_vlm_model=cloud_vlm_model,
        )
        self._llm = LlamaLLMClient(
            base_url=settings.NIM_LLM_URL,
            mock_enabled=mock_enabled,
            anthropic_api_key=getattr(settings, "ANTHROPIC_API_KEY", None),
            nvidia_api_key=resolved_nvidia_key,
            cloud_url=cloud_url,
            cloud_llm_model=cloud_llm_model,
            local_llm_model=getattr(settings, "NIM_LOCAL_LLM_MODEL", "meta/llama3-70b-instruct"),
        )

        logger.info("NIM Service Manager initialized")
        logger.info(f"  VISTA-3D : {settings.NIM_VISTA3D_URL}")
        logger.info(f"  MAISI    : {settings.NIM_MAISI_URL}")
        logger.info(f"  VILA-M3  : {settings.NIM_VILAM3_URL}")
        logger.info(f"  LLM      : {settings.NIM_LLM_URL}")
        logger.info(f"  Cloud URL: {cloud_url}")
        logger.info(f"  Cloud LLM: {cloud_llm_model}")
        logger.info(f"  Cloud VLM: {cloud_vlm_model}")
        logger.info(f"  NVIDIA API key: {'set' if resolved_nvidia_key else 'not set'}")
        logger.info(f"  Mock fallback: {mock_enabled}")

    # ── Typed properties ──

    @property
    def vista3d(self) -> VISTA3DClient:
        """VISTA-3D segmentation client."""
        return self._vista3d

    @property
    def maisi(self) -> MAISIClient:
        """MAISI synthetic CT generation client."""
        return self._maisi

    @property
    def vilam3(self) -> VILAM3Client:
        """VILA-M3 visual language model client."""
        return self._vilam3

    @property
    def llm(self) -> LlamaLLMClient:
        """LLM (Llama-3 / Cloud NIM / Claude) client."""
        return self._llm

    # ── Service management ──

    def check_all_services(self) -> Dict[str, str]:
        """Check health of all NIM services.

        Returns:
            Dictionary mapping service name to status string.
            Status is one of: "available", "cloud", "anthropic", "mock", "unavailable".
        """
        status = {
            "vista3d": self._vista3d.get_status(),
            "maisi": self._maisi.get_status(),
            "vila_m3": self._vilam3.get_status(),
            "llm": self._llm.get_status(),
        }

        available_count = sum(1 for s in status.values() if s == "available")
        cloud_count = sum(1 for s in status.values() if s == "cloud")
        mock_count = sum(1 for s in status.values() if s == "mock")
        logger.info(
            f"NIM service status: {available_count} available, "
            f"{cloud_count} cloud, "
            f"{mock_count} mock, "
            f"{len(status) - available_count - cloud_count - mock_count} other"
        )

        return status

    def get_available_services(self) -> List[str]:
        """Get names of services that are live (not mock or unavailable).

        Returns:
            List of service names with status "available" or "cloud".
        """
        status = self.check_all_services()
        return [name for name, state in status.items() if state in ("available", "cloud")]

    def get_client(self, service_name: str) -> BaseNIMClient:
        """Get a NIM client by service name.

        Args:
            service_name: One of "vista3d", "maisi", "vila_m3", "llm".

        Returns:
            The corresponding BaseNIMClient subclass instance.

        Raises:
            KeyError: If service_name is not recognized.
        """
        clients = {
            "vista3d": self._vista3d,
            "maisi": self._maisi,
            "vila_m3": self._vilam3,
            "llm": self._llm,
        }

        if service_name not in clients:
            raise KeyError(
                f"Unknown NIM service: {service_name!r}. "
                f"Available: {list(clients.keys())}"
            )

        return clients[service_name]
