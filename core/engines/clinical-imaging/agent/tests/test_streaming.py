"""Tests for NVIDIA Holoscan SDK real-time streaming integration.

All tests work without Holoscan SDK or OpenCV installed.
Uses numpy arrays as mock frames for reproducible testing.

Author: Adam Jones
Date: April 2026
"""

import numpy as np
import pytest

from src.streaming.holoscan_manager import (
    HoloscanStreamManager,
    StreamConfig,
    StreamStatus,
)
from src.streaming.base_streaming_workflow import BaseStreamingWorkflow
from src.streaming.ultrasound_realtime import UltrasoundRealtimeWorkflow
from src.streaming.endoscopy_realtime import EndoscopyRealtimeWorkflow


# ===================================================================
# HELPERS
# ===================================================================


def _make_frame(height: int = 480, width: int = 640) -> np.ndarray:
    """Generate a random mock video frame."""
    return np.random.randint(0, 255, (height, width, 3), dtype=np.uint8)


# ===================================================================
# HOLOSCAN MANAGER — GRACEFUL FALLBACK
# ===================================================================


class TestHoloscanGracefulFallback:
    """Verify the manager initializes without Holoscan SDK."""

    def test_holoscan_graceful_fallback(self):
        """HoloscanStreamManager works when holoscan is not installed."""
        manager = HoloscanStreamManager()
        # Should not raise — gracefully detects missing SDK
        assert isinstance(manager, HoloscanStreamManager)
        # holoscan_available is a bool (True if installed, False otherwise)
        assert isinstance(manager.holoscan_available, bool)


# ===================================================================
# STREAM CONFIG & STATUS DEFAULTS
# ===================================================================


class TestStreamConfig:
    """Verify StreamConfig default values."""

    def test_stream_config_defaults(self):
        """StreamConfig has correct default values."""
        config = StreamConfig()
        assert config.source_type == "file"
        assert config.target_fps == 30
        assert config.display is False
        assert config.record is False
        assert config.output_dir == "data/streaming_output"

    def test_stream_config_custom(self):
        """StreamConfig accepts custom values."""
        config = StreamConfig(
            source_type="rtsp",
            source_path="rtsp://192.168.1.100:554/live",
            target_fps=15,
            display=True,
        )
        assert config.source_type == "rtsp"
        assert config.source_path == "rtsp://192.168.1.100:554/live"
        assert config.target_fps == 15


class TestStreamStatus:
    """Verify StreamStatus default values."""

    def test_stream_status_defaults(self):
        """StreamStatus initializes with correct defaults."""
        status = StreamStatus()
        assert status.active is False
        assert status.frames_processed == 0
        assert status.fps_actual == 0.0
        assert status.findings_count == 0
        assert status.elapsed_seconds == 0.0
        assert status.last_finding is None


# ===================================================================
# ULTRASOUND REALTIME WORKFLOW
# ===================================================================


class TestUltrasoundRealtime:
    """Tests for UltrasoundRealtimeWorkflow mock mode."""

    def test_ultrasound_mock_frame(self):
        """Process a single frame in mock mode and get a valid result."""
        workflow = UltrasoundRealtimeWorkflow(mock_mode=True)
        frame = _make_frame()
        result = workflow.process_frame(frame, frame_index=0)

        assert isinstance(result, dict)
        assert "detections" in result
        assert "frame_index" in result
        assert result["frame_index"] == 0
        assert isinstance(result["detections"], list)

    def test_ultrasound_detection_frequency(self):
        """Detections appear periodically, not every frame."""
        workflow = UltrasoundRealtimeWorkflow(
            mock_mode=True, detection_interval=90
        )
        frame = _make_frame()

        detection_frames = []
        empty_frames = []

        for i in range(300):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection_frames.append(i)
            else:
                empty_frames.append(i)

        # Should have some detections but many more empty frames
        assert len(detection_frames) > 0, "Expected at least one detection"
        assert len(empty_frames) > len(detection_frames), (
            "Most frames should be detection-free"
        )
        # Frame 0 should not have a detection (interval check starts > 0)
        assert 0 not in detection_frames

    def test_ultrasound_detection_has_bbox(self):
        """Detections include bounding box coordinates."""
        workflow = UltrasoundRealtimeWorkflow(
            mock_mode=True, detection_interval=90
        )
        frame = _make_frame()

        # Find a frame with a detection
        detection = None
        for i in range(300):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection = result["detections"][0]
                break

        assert detection is not None, "Expected at least one detection in 300 frames"
        assert "bbox" in detection
        bbox = detection["bbox"]
        assert "x_min" in bbox
        assert "y_min" in bbox
        assert "x_max" in bbox
        assert "y_max" in bbox
        # Coordinates should be in [0, 1] range (normalized)
        assert 0.0 <= bbox["x_min"] <= 1.0
        assert 0.0 <= bbox["y_min"] <= 1.0
        assert bbox["x_max"] > bbox["x_min"]
        assert bbox["y_max"] > bbox["y_min"]

    def test_ultrasound_detection_has_scoring(self):
        """Detections include BI-RADS or TI-RADS scoring."""
        workflow = UltrasoundRealtimeWorkflow(
            mock_mode=True, scoring_system="birads"
        )
        frame = _make_frame()

        detection = None
        for i in range(300):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection = result["detections"][0]
                break

        assert detection is not None
        assert "score" in detection
        assert "scoring_system" in detection
        assert detection["scoring_system"] == "birads"
        assert "confidence" in detection
        assert 0.0 <= detection["confidence"] <= 1.0

    def test_ultrasound_workflow_info(self):
        """Workflow returns correct metadata."""
        workflow = UltrasoundRealtimeWorkflow(mock_mode=True)
        info = workflow.get_workflow_info()

        assert info["name"] == "ultrasound_realtime"
        assert info["modality"] == "ultrasound"
        assert info["target_fps"] == 30
        assert info["type"] == "streaming"
        assert info["mock_mode"] is True

    def test_ultrasound_tirads_mode(self):
        """TI-RADS scoring system produces thyroid-specific detections."""
        workflow = UltrasoundRealtimeWorkflow(
            mock_mode=True, scoring_system="tirads"
        )
        frame = _make_frame()

        detection = None
        for i in range(300):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection = result["detections"][0]
                break

        assert detection is not None
        assert detection["scoring_system"] == "tirads"
        assert detection["label"] == "thyroid_nodule"
        assert detection["score"].startswith("TR")


# ===================================================================
# ENDOSCOPY REALTIME WORKFLOW
# ===================================================================


class TestEndoscopyRealtime:
    """Tests for EndoscopyRealtimeWorkflow mock mode."""

    def test_endoscopy_mock_frame(self):
        """Process a single frame in mock mode and get a valid result."""
        workflow = EndoscopyRealtimeWorkflow(mock_mode=True)
        frame = _make_frame()
        result = workflow.process_frame(frame, frame_index=0)

        assert isinstance(result, dict)
        assert "detections" in result
        assert "frame_index" in result
        assert "alert" in result
        assert result["frame_index"] == 0

    def test_endoscopy_polyp_detection(self):
        """Polyp detections include Paris classification."""
        workflow = EndoscopyRealtimeWorkflow(
            mock_mode=True, detection_interval=120
        )
        frame = _make_frame()

        detection = None
        for i in range(500):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection = result["detections"][0]
                break

        assert detection is not None, "Expected at least one detection in 500 frames"
        assert "paris_class" in detection
        assert detection["paris_class"].startswith("0-")
        assert "paris_description" in detection
        assert "label" in detection
        assert detection["label"] == "polyp"

    def test_endoscopy_size_estimation(self):
        """Polyp detections include size_mm estimate."""
        workflow = EndoscopyRealtimeWorkflow(
            mock_mode=True, detection_interval=120
        )
        frame = _make_frame()

        detection = None
        for i in range(500):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection = result["detections"][0]
                break

        assert detection is not None
        assert "size_mm" in detection
        assert isinstance(detection["size_mm"], float)
        assert detection["size_mm"] > 0.0

    def test_endoscopy_nice_classification(self):
        """Polyp detections include NICE classification."""
        workflow = EndoscopyRealtimeWorkflow(
            mock_mode=True, detection_interval=120
        )
        frame = _make_frame()

        detection = None
        for i in range(500):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection = result["detections"][0]
                break

        assert detection is not None
        assert "nice_class" in detection
        assert detection["nice_class"] in ("Type 1", "Type 2", "Type 3")

    def test_endoscopy_location(self):
        """Polyp detections include anatomical location."""
        workflow = EndoscopyRealtimeWorkflow(
            mock_mode=True, detection_interval=120
        )
        frame = _make_frame()

        detection = None
        for i in range(500):
            result = workflow.process_frame(frame, frame_index=i)
            if result["detections"]:
                detection = result["detections"][0]
                break

        assert detection is not None
        assert "location" in detection
        assert detection["location"] in [
            "cecum", "ascending_colon", "hepatic_flexure",
            "transverse_colon", "splenic_flexure", "descending_colon",
            "sigmoid_colon", "rectum",
        ]


# ===================================================================
# HOLOSCAN STREAM MANAGER
# ===================================================================


class TestStreamManager:
    """Tests for HoloscanStreamManager operations."""

    def test_process_single_frame(self):
        """HoloscanStreamManager.process_single_frame works with mock workflow."""
        manager = HoloscanStreamManager()
        workflow = UltrasoundRealtimeWorkflow(mock_mode=True)
        frame = _make_frame()

        result = manager.process_single_frame(frame, workflow)

        assert isinstance(result, dict)
        assert "detections" in result
        assert "frame_index" in result
        assert result["frame_index"] == 0

    def test_list_streams_empty(self):
        """No active streams initially."""
        manager = HoloscanStreamManager()
        streams = manager.list_streams()

        assert isinstance(streams, list)
        assert len(streams) == 0

    def test_stream_manager_status(self):
        """get_status returns inactive status for unknown stream."""
        manager = HoloscanStreamManager()
        status = manager.get_status("nonexistent-stream")

        assert isinstance(status, StreamStatus)
        assert status.active is False
        assert status.frames_processed == 0
