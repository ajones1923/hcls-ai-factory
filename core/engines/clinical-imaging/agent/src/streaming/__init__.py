"""NVIDIA Holoscan SDK integration for real-time streaming inference.

Provides HoloscanStreamManager for managing real-time video analysis
pipelines, plus concrete workflows for ultrasound and endoscopy.

Apache 2.0 Licensed. ARM64 native (designed for Jetson -> IGX -> DGX).
"""

from src.streaming.holoscan_manager import (
    HoloscanStreamManager,
    StreamConfig,
    StreamStatus,
)
from src.streaming.base_streaming_workflow import BaseStreamingWorkflow
from src.streaming.ultrasound_realtime import UltrasoundRealtimeWorkflow
from src.streaming.endoscopy_realtime import EndoscopyRealtimeWorkflow

__all__ = [
    "HoloscanStreamManager",
    "StreamConfig",
    "StreamStatus",
    "BaseStreamingWorkflow",
    "UltrasoundRealtimeWorkflow",
    "EndoscopyRealtimeWorkflow",
]
