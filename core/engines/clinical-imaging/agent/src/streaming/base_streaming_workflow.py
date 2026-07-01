"""Base class for real-time streaming analysis workflows.

All streaming workflows process individual video frames at target FPS
(typically 30fps = 33ms per frame budget). This differs from the static
BaseImagingWorkflow which processes whole volumes/images.

Apache 2.0 Licensed. ARM64 native (designed for Jetson -> IGX -> DGX).
"""

from abc import ABC, abstractmethod
from typing import Dict, List

import numpy as np
from loguru import logger


class BaseStreamingWorkflow(ABC):
    """Base class for real-time streaming analysis workflows.

    Subclasses implement process_frame() to analyze a single video frame
    within the per-frame time budget (e.g., <33ms for 30fps).
    """

    WORKFLOW_NAME: str = "base_streaming"
    TARGET_FPS: int = 30
    MODALITY: str = ""
    MODELS_USED: List[str] = []

    def __init__(self, mock_mode: bool = True):
        self.mock_mode = mock_mode
        logger.info(
            f"Initialized {self.WORKFLOW_NAME} streaming workflow "
            f"(mock={mock_mode}, target_fps={self.TARGET_FPS})"
        )

    @abstractmethod
    def process_frame(self, frame: np.ndarray, frame_index: int) -> Dict:
        """Process a single video frame. Must complete in <33ms for 30fps.

        Args:
            frame: BGR image array of shape (H, W, 3), dtype uint8.
            frame_index: Zero-based frame counter.

        Returns:
            Dict with at minimum:
                - detections: list of detection dicts (may be empty)
                - frame_index: echo of the input frame_index
        """
        ...

    @abstractmethod
    def _mock_frame_result(self, frame_index: int) -> Dict:
        """Return a mock result for a single frame.

        Used when mock_mode=True or when real inference is unavailable.
        Should produce clinically realistic detection patterns.
        """
        ...

    def get_workflow_info(self) -> Dict:
        """Return metadata about this streaming workflow."""
        return {
            "name": self.WORKFLOW_NAME,
            "modality": self.MODALITY,
            "target_fps": self.TARGET_FPS,
            "models_used": self.MODELS_USED,
            "mock_mode": self.mock_mode,
            "type": "streaming",
        }
