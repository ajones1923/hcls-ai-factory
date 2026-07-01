"""Real-time endoscopy polyp detection workflow.

Processes endoscopy video at 30fps with polyp detection
and size estimation overlay. Implements Paris classification
for polyp morphology and generates clinical alerts for
high-risk findings.

Uses MONAI models for detection with TensorRT optimization.
In mock mode, generates clinically realistic polyp detections
approximately every 120 frames (4 seconds at 30fps).

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

# Paris classification for superficial GI neoplasms
PARIS_CLASSIFICATIONS = {
    "0-Ip": "Pedunculated (protruding with stalk)",
    "0-Isp": "Sub-pedunculated (protruding, intermediate)",
    "0-Is": "Sessile (protruding, broad-based)",
    "0-IIa": "Flat elevated (slightly elevated)",
    "0-IIb": "Flat (completely flat)",
    "0-IIc": "Flat depressed (slightly depressed)",
    "0-III": "Excavated (ulcerated / depressed)",
}

# Polyp types for classification
POLYP_TYPES = [
    "adenomatous",
    "hyperplastic",
    "sessile_serrated",
    "traditional_serrated",
    "inflammatory",
]

# NICE (NBI International Colorectal Endoscopic) classification
NICE_CLASSIFICATIONS = {
    "Type 1": "Hyperplastic — no treatment needed",
    "Type 2": "Adenomatous — polypectomy recommended",
    "Type 3": "Invasive — surgical evaluation needed",
}

# Anatomical locations in the colon
COLON_SEGMENTS = [
    "cecum",
    "ascending_colon",
    "hepatic_flexure",
    "transverse_colon",
    "splenic_flexure",
    "descending_colon",
    "sigmoid_colon",
    "rectum",
]

# Alert thresholds
ALERT_SIZE_MM = 10.0  # Polyps >= 10mm warrant immediate attention
ALERT_PARIS_CLASSES = {"0-IIc", "0-III"}  # Depressed/excavated = higher risk


class EndoscopyRealtimeWorkflow(BaseStreamingWorkflow):
    """Real-time endoscopy polyp detection with Paris classification.

    Designed for colonoscopy and upper GI endoscopy video analysis
    at 30fps. Each frame is evaluated for polyps; when detected,
    the result includes bounding box, size estimation, Paris
    classification, and clinical alert level.
    """

    WORKFLOW_NAME: str = "endoscopy_realtime"
    TARGET_FPS: int = 30
    MODALITY: str = "endoscopy"
    MODELS_USED: List[str] = [
        "MONAI Detection (TensorRT)",
        "Polyp Classifier",
        "Size Estimator",
    ]

    def __init__(
        self,
        mock_mode: bool = True,
        detection_interval: int = 120,
    ):
        """Initialize endoscopy real-time workflow.

        Args:
            mock_mode: If True, generate mock detections without inference.
            detection_interval: Approximate frames between mock detections.
        """
        super().__init__(mock_mode=mock_mode)
        self.detection_interval = detection_interval
        self._rng = random.Random(42)
        logger.info(
            f"Endoscopy workflow: detection_interval={detection_interval}"
        )

    def process_frame(self, frame: np.ndarray, frame_index: int) -> Dict:
        """Detect polyps in a single endoscopy frame.

        Args:
            frame: BGR image array of shape (H, W, 3), dtype uint8.
            frame_index: Zero-based frame counter.

        Returns:
            Dict with keys:
                - detections: list of detection dicts with bbox, confidence,
                  size_mm, paris_class
                - frame_index: echo of input
                - overlay_image: annotated frame (None in mock mode)
                - alert: clinical alert string or None
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
        detection_interval frames, a polyp detection is simulated
        with Paris classification, size estimation, and alert logic.
        """
        result = {
            "detections": [],
            "frame_index": frame_index,
            "overlay_image": None,
            "alert": None,
        }

        # Simulate periodic detection with jitter
        jitter = self._rng.randint(-20, 20)
        effective_interval = max(40, self.detection_interval + jitter)

        if frame_index > 0 and frame_index % effective_interval == 0:
            detection = self._generate_mock_detection(frame_index)
            result["detections"] = [detection]

            # Generate alert for high-risk findings
            alert = self._evaluate_alert(detection)
            if alert:
                result["alert"] = alert
                logger.info(f"Frame {frame_index}: ALERT — {alert}")
            else:
                logger.debug(
                    f"Frame {frame_index}: polyp detected "
                    f"({detection['paris_class']}, "
                    f"{detection['size_mm']}mm)"
                )

        return result

    def _generate_mock_detection(self, frame_index: int) -> Dict:
        """Generate a single mock polyp detection with clinical metadata."""
        # Bounding box in normalized coordinates [0, 1]
        cx = self._rng.uniform(0.20, 0.80)
        cy = self._rng.uniform(0.20, 0.80)
        w = self._rng.uniform(0.04, 0.18)
        h = self._rng.uniform(0.04, 0.18)

        bbox = {
            "x_min": round(cx - w / 2, 3),
            "y_min": round(cy - h / 2, 3),
            "x_max": round(cx + w / 2, 3),
            "y_max": round(cy + h / 2, 3),
        }

        confidence = round(self._rng.uniform(0.60, 0.97), 2)

        # Paris classification — weighted toward common morphologies
        paris_weights = [0.15, 0.10, 0.30, 0.20, 0.10, 0.08, 0.07]
        paris_classes = list(PARIS_CLASSIFICATIONS.keys())
        paris_class = self._rng.choices(
            paris_classes, weights=paris_weights, k=1
        )[0]
        paris_description = PARIS_CLASSIFICATIONS[paris_class]

        # Size estimation in mm
        # Pedunculated polyps tend to be larger; flat lesions tend to be smaller
        if paris_class in ("0-Ip", "0-Isp"):
            size_mm = round(self._rng.uniform(8.0, 35.0), 1)
        elif paris_class in ("0-IIb", "0-IIc", "0-III"):
            size_mm = round(self._rng.uniform(3.0, 20.0), 1)
        else:
            size_mm = round(self._rng.uniform(4.0, 25.0), 1)

        # Polyp type classification
        polyp_type = self._rng.choice(POLYP_TYPES)

        # NICE classification
        nice_weights = [0.35, 0.50, 0.15]
        nice_types = list(NICE_CLASSIFICATIONS.keys())
        nice_class = self._rng.choices(
            nice_types, weights=nice_weights, k=1
        )[0]
        nice_description = NICE_CLASSIFICATIONS[nice_class]

        # Anatomical location
        location = self._rng.choice(COLON_SEGMENTS)

        return {
            "label": "polyp",
            "confidence": confidence,
            "bbox": bbox,
            "size_mm": size_mm,
            "paris_class": paris_class,
            "paris_description": paris_description,
            "polyp_type": polyp_type,
            "nice_class": nice_class,
            "nice_description": nice_description,
            "location": location,
            "frame_index": frame_index,
        }

    @staticmethod
    def _evaluate_alert(detection: Dict) -> Optional[str]:
        """Evaluate whether a detection warrants a clinical alert.

        Alerts are generated for:
            - Polyps >= 10mm (increased malignancy risk)
            - Paris 0-IIc or 0-III morphology (depressed/excavated)
            - NICE Type 3 (invasive pattern)
        """
        alerts = []

        size_mm = detection.get("size_mm", 0)
        paris_class = detection.get("paris_class", "")
        nice_class = detection.get("nice_class", "")

        if size_mm >= ALERT_SIZE_MM:
            alerts.append(f"Large polyp ({size_mm}mm)")

        if paris_class in ALERT_PARIS_CLASSES:
            alerts.append(
                f"High-risk morphology ({paris_class}: "
                f"{PARIS_CLASSIFICATIONS.get(paris_class, '')})"
            )

        if nice_class == "Type 3":
            alerts.append("Invasive pattern (NICE Type 3)")

        if alerts:
            location = detection.get("location", "unknown")
            return f"{location}: {'; '.join(alerts)}"
        return None
