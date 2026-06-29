"""Base NIM client with health check, retry, and mock fallback."""

import time
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

import requests
from loguru import logger
from tenacity import retry, stop_after_attempt, wait_exponential


class BaseNIMClient(ABC):
    """Abstract base for all NIM service clients."""

    def __init__(self, base_url: str, service_name: str, mock_enabled: bool = True):
        self.base_url = base_url.rstrip("/")
        self.service_name = service_name
        self.mock_enabled = mock_enabled
        self._available: Optional[bool] = None
        self._last_check: float = 0
        self._check_interval = 30.0  # seconds

    def health_check(self) -> bool:
        """Ping the NIM health endpoint."""
        try:
            resp = requests.get(f"{self.base_url}/v1/health/ready", timeout=5)
            return resp.status_code == 200
        except Exception:
            return False

    def is_available(self) -> bool:
        """Cached availability check."""
        now = time.time()
        if self._available is None or (now - self._last_check) > self._check_interval:
            self._available = self.health_check()
            self._last_check = now
            if self._available:
                logger.info(f"NIM {self.service_name} is available at {self.base_url}")
            else:
                logger.warning(f"NIM {self.service_name} unavailable at {self.base_url}")
        return self._available

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=10))
    def _request(self, endpoint: str, payload: Dict, timeout: int = 120) -> Dict:
        """HTTP POST with retry logic."""
        url = f"{self.base_url}{endpoint}"
        resp = requests.post(url, json=payload, timeout=timeout)
        resp.raise_for_status()
        return resp.json()

    @abstractmethod
    def _mock_response(self, **kwargs) -> Any:
        """Return mock response for demo mode."""
        ...

    def _invoke_or_mock(
        self, endpoint: str, payload: Dict, timeout: int = 120, **mock_kwargs
    ) -> Any:
        """Try real NIM, fall back to mock if unavailable and mock_enabled."""
        if self.is_available():
            try:
                return self._request(endpoint, payload, timeout)
            except Exception as e:
                logger.error(f"NIM {self.service_name} request failed: {e}")
                if self.mock_enabled:
                    logger.warning(f"Falling back to mock for {self.service_name}")
                    return self._mock_response(**mock_kwargs)
                raise
        elif self.mock_enabled:
            logger.info(f"Using mock for {self.service_name} (service unavailable)")
            return self._mock_response(**mock_kwargs)
        else:
            raise ConnectionError(
                f"NIM {self.service_name} unavailable and mock disabled"
            )

    def get_status(self) -> str:
        """Return service status string."""
        if self.is_available():
            return "available"
        elif self.mock_enabled:
            return "mock"
        return "unavailable"
