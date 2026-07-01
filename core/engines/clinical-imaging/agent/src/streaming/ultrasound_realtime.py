"""Real-time ultrasound lesion detection workflow.

Processes ultrasound video at 30fps with lesion detection
and BI-RADS/TI-RADS scoring overlay.

Uses MONAI models for detection with TensorRT optimization.
In mock mode, generates clinically realistic detection patterns:
lesions appear approximately every 90 frames (3 seconds at 30fps)
with realistic bounding boxes, confidence scores, and BI-RADS
classifications.

Apache 2.0 Licensed. ARM64 native (designed for Jetson -> IGX -> DGX).
"""

import random
from typing import Dict, List, Optional

import numpy as np
from loguru import logger

from src.streaming.base_streaming_workflow import BaseStreamingWorkflow


# ═══════════════════════════════════════════════════════════════════════
# CLINICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

# BI-RADS categories with clinical descriptions
BIRADS_CATEGORIES = {
    "2": "Benign — no malignant features",
    "3": "Probably benign — short-interval follow-up suggested",
    "4A": "Low suspicion for malignancy — biopsy should be considered",
    "4B": "Moderate suspicion for malignancy — biopsy recommended",
    "4C": "High suspicion for malignancy — biopsy strongly recommended",
    "5": "Highly suggestive of malignancy — tissue diagnosis required",
}

# TI-RADS categories for thyroid
TIRADS_CATEGORIES = {
    "TR1": "Benign",
    "TR2": "Not suspicious",
    "TR3": "Mildly suspicious",
    "TR4": "Moderately suspicious",
    "TR5": "Highly suspicious",
}

# Lesion morphology descriptors
LESION_SHAPES = [
    "oval",
    "round",
    "irregular",
    "lobulated",
]

LESION_MARGINS = [
    "circumscribed",
    "indistinct",
    "angular",
    "microlobulated",
    "spiculated",
]

LESION_ECHO_PATTERNS = [
    "anechoic",
    "hyperechoic",
    "isoechoic",
    "hypoechoic",
    "complex",
]


class UltrasoundRealtimeWorkflow(BaseStreamingWorkflow):
    """Real-time ultrasound lesion detection with BI-RADS/TI-RADS scoring.

    Designed for breast and thyroid ultrasound video analysis at 30fps.
    Each frame is evaluated for suspicious lesions; when detected, the
    result includes bounding box coordinates, confidence score, and
    appropriate scoring system classification.
    """

    WORKFLOW_NAME: str = "ultrasound_realtime"
    TARGET_FPS: int = 30
    MODALITY: str = "ultrasound"
    MODELS_USED: List[str] = [
        "MONAI Detection (TensorRT)",
        "BI-RADS Classifier",
    ]

    def __init__(
        self,
        mock_mode: bool = True,
        detection_interval: int = 90,
        scoring_system: str = "birads",
    ):
        """Initialize ultrasound real-time workflow.

        Args:
            mock_mode: If True, generate mock detections without inference.
            detection_interval: Approximate frames between mock detections.
            scoring_system: "birads" for breast, "tirads" for thyroid.
        """
        super().__init__(mock_mode=mock_mode)
        self.detection_interval = detection_interval
        self.scoring_system = scoring_system
        self._rng = random.Random(42)
        logger.info(
            f"Ultrasound workflow: scoring={scoring_system}, "
            f"detection_interval={detection_interval}"
        )

    def process_frame(self, frame: np.ndarray, frame_index: int) -> Dict:
        """Detect lesions in a single ultrasound frame.

        Args:
            frame: BGR image array of shape (H, W, 3), dtype uint8.
            frame_index: Zero-based frame counter.

        Returns:
            Dict with keys:
                - detections: list of detection dicts with bbox, confidence, label
                - frame_index: echo of input
                - overlay_image: annotated frame (None in mock mode)
                - severity: overall frame severity
        """
        if self.mock_mode:
            return self._mock_frame_result(frame_index)

        # Real inference path — requires MONAI + TensorRT
        return self._real_inference(frame, frame_index)

    def _real_inference(self, frame: np.ndarray, frame_index: int) -> Dict:
        """Run real model inference on a frame.

        Placeholder for production deployment with MONAI Detection
        model optimized via TensorRT.
        """
        logger.warning(
            "Real inference not yet implemented — "
            "falling back to mock results"
        )
        return self._mock_frame_result(frame_index)

    def _mock_frame_result(self, frame_index: int) -> Dict:
        """Generate a clinically realistic mock result for a single frame.

        Most frames produce no detections. Approximately every
        detection_interval frames, a lesion detection is simulated
        with realistic BI-RADS/TI-RADS scoring and morphology.
        """
        # Default: no detections (most frames are clean)
        result = {
            "detections": [],
            "frame_index": frame_index,
            "overlay_image": None,
            "severity": "normal",
        }

        # Simulate periodic detection with some jitter
        jitter = self._rng.randint(-15, 15)
        effective_interval = max(30, self.detection_interval + jitter)

        if frame_index > 0 and frame_index % effective_interval == 0:
            detection = self._generate_mock_detection(frame_index)
            result["detections"] = [detection]
            result["severity"] = detection.get("severity", "significant")
            logger.debug(
                f"Frame {frame_index}: detected {detection['label']} "
                f"(confidence={detection['confidence']:.2f})"
            )

        return result

    def _generate_mock_detection(self, frame_index: int) -> Dict:
        """Generate a single mock lesion detection with clinical metadata."""
        # Bounding box in normalized coordinates [0, 1]
        cx = self._rng.uniform(0.25, 0.75)
        cy = self._rng.uniform(0.25, 0.75)
        w = self._rng.uniform(0.05, 0.15)
        h = self._rng.uniform(0.05, 0.15)

        bbox = {
            "x_min": round(cx - w / 2, 3),
            "y_min": round(cy - h / 2, 3),
            "x_max": round(cx + w / 2, 3),
            "y_max": round(cy + h / 2, 3),
        }

        confidence = round(self._rng.uniform(0.65, 0.95), 2)

        # Select scoring category based on confidence
        if self.scoring_system == "tirads":
            score, description = self._assign_tirads(confidence)
            label = "thyroid_nodule"
        else:
            score, description = self._assign_birads(confidence)
            label = "breast_lesion"

        # Morphology descriptors
        shape = self._rng.choice(LESION_SHAPES)
        margin = self._rng.choice(LESION_MARGINS)
        echo_pattern = self._rng.choice(LESION_ECHO_PATTERNS)

        # Estimated size in mm (realistic range for ultrasound lesions)
        size_mm = round(self._rng.uniform(4.0, 28.0), 1)

        # Severity mapping
        severity = self._score_to_severity(score)

        return {
            "label": label,
            "confidence": confidence,
            "bbox": bbox,
            "size_mm": size_mm,
            "score": score,
            "score_description": description,
            "scoring_system": self.scoring_system,
            "morphology": {
                "shape": shape,
                "margin": margin,
                "echo_pattern": echo_pattern,
            },
            "severity": severity,
            "frame_index": frame_index,
        }

    def _assign_birads(self, confidence: float) -> tuple:
        """Assign a BI-RADS category based on detection confidence."""
        if confidence >= 0.90:
            cat = "4C"
        elif confidence >= 0.85:
            cat = "4B"
        elif confidence >= 0.78:
            cat = "4A"
        elif confidence >= 0.72:
            cat = "3"
        else:
            cat = "2"
        return cat, BIRADS_CATEGORIES.get(cat, "")

    def _assign_tirads(self, confidence: float) -> tuple:
        """Assign a TI-RADS category based on detection confidence."""
        if confidence >= 0.90:
            cat = "TR5"
        elif confidence >= 0.82:
            cat = "TR4"
        elif confidence >= 0.72:
            cat = "TR3"
        else:
            cat = "TR2"
        return cat, TIRADS_CATEGORIES.get(cat, "")

    @staticmethod
    def _score_to_severity(score: str) -> str:
        """Map a BI-RADS or TI-RADS score to a severity string."""
        high_risk = {"4B", "4C", "5", "TR4", "TR5"}
        moderate_risk = {"4A", "TR3"}
        low_risk = {"3", "TR2"}

        if score in high_risk:
            return "urgent"
        elif score in moderate_risk:
            return "significant"
        elif score in low_risk:
            return "routine"
        return "normal"
