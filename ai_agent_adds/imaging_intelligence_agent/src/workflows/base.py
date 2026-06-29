"""Base imaging workflow with mock support."""

import time
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import FindingSeverity, WorkflowResult, WorkflowStatus


class BaseImagingWorkflow(ABC):
    """Abstract base class for all imaging inference workflows."""

    WORKFLOW_NAME: str = "base"
    TARGET_LATENCY_SEC: float = 60.0
    MODALITY: str = ""
    BODY_REGION: str = ""
    MODELS_USED: List[str] = []

    def __init__(self, mock_mode: bool = True, nim_clients: Optional[Dict] = None):
        self.mock_mode = mock_mode
        self.nim_clients = nim_clients or {}
        logger.info(f"Initialized {self.WORKFLOW_NAME} workflow (mock={mock_mode})")

    @abstractmethod
    def preprocess(self, input_path: str) -> Any:
        """Load and preprocess input image/volume."""
        ...

    @abstractmethod
    def infer(self, preprocessed: Any) -> Dict:
        """Run model inference (real or mock)."""
        ...

    @abstractmethod
    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Extract measurements, classifications, findings from inference."""
        ...

    @abstractmethod
    def _mock_inference(self) -> Dict:
        """Return realistic mock inference result."""
        ...

    def run(self, input_path: str = "") -> WorkflowResult:
        """Orchestrate the full workflow pipeline."""
        start = time.time()
        try:
            if self.mock_mode:
                logger.info(f"Running {self.WORKFLOW_NAME} in mock mode")
                raw = self._mock_inference()
            else:
                preprocessed = self.preprocess(input_path)
                raw = self.infer(preprocessed)

            result = self.postprocess(raw)
            result.inference_time_ms = (time.time() - start) * 1000
            result.is_mock = self.mock_mode
            result.workflow_name = self.WORKFLOW_NAME
            return result
        except Exception as e:
            logger.error(f"Workflow {self.WORKFLOW_NAME} failed: {e}")
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.FAILED,
                inference_time_ms=(time.time() - start) * 1000,
                is_mock=self.mock_mode,
            )

    def get_workflow_info(self) -> Dict:
        """Return metadata about this workflow."""
        return {
            "name": self.WORKFLOW_NAME,
            "modality": self.MODALITY,
            "body_region": self.BODY_REGION,
            "target_latency_sec": self.TARGET_LATENCY_SEC,
            "models_used": self.MODELS_USED,
            "mock_mode": self.mock_mode,
        }
