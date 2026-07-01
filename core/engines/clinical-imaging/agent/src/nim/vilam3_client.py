"""VILA-M3 NIM client for medical visual language understanding.

VILA-M3 (Visual Language Model for Medical 3D) provides visual question
answering, report generation, and zero-shot classification on medical
images using a Llama-3 backbone with visual encoding.

Fallback chain: local VILA-M3 NIM -> NVIDIA cloud Llama 3.2 Vision -> mock
"""

import base64
import random
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.models import VLMResponse

from .base import BaseNIMClient


class VILAM3Client(BaseNIMClient):
    """Client for NVIDIA VILA-M3 NIM visual language model service.

    VILA-M3 uses an OpenAI-compatible API format with a Llama-3
    backbone for multimodal medical image understanding. Capabilities:
      - Visual question answering on radiology images
      - Structured radiology report generation
      - Zero-shot image classification

    Fallback chain:
      1. Local VILA-M3 NIM endpoint
      2. NVIDIA Cloud NIM (Llama-3.2-11B-Vision via integrate.api.nvidia.com)
      3. Mock response (if mock_enabled)
    """

    def __init__(
        self,
        base_url: str,
        mock_enabled: bool = True,
        nvidia_api_key: Optional[str] = None,
        cloud_url: str = "https://integrate.api.nvidia.com/v1",
        cloud_vlm_model: str = "meta/llama-3.2-11b-vision-instruct",
    ):
        super().__init__(base_url, service_name="vila_m3", mock_enabled=mock_enabled)
        self.nvidia_api_key = nvidia_api_key
        self.cloud_url = cloud_url.rstrip("/")
        self.cloud_vlm_model = cloud_vlm_model
        self._cloud_available: Optional[bool] = None

    def cloud_health_check(self) -> bool:
        """Check NVIDIA cloud VLM NIM availability."""
        if not self.nvidia_api_key:
            return False
        try:
            resp = requests.get(
                f"{self.cloud_url}/models",
                headers={"Authorization": f"Bearer {self.nvidia_api_key}"},
                timeout=10,
            )
            return resp.status_code == 200
        except Exception:
            return False

    def is_cloud_available(self) -> bool:
        """Check if NVIDIA cloud VLM NIM is available (cached)."""
        import time
        now = time.time()
        if self._cloud_available is None or (now - self._last_check) > self._check_interval:
            self._cloud_available = self.cloud_health_check()
            if self._cloud_available:
                logger.info(f"NVIDIA Cloud NIM VLM available at {self.cloud_url}")
            else:
                logger.debug("NVIDIA Cloud NIM VLM not available")
        return self._cloud_available

    def _encode_image(self, image_path: str) -> str:
        """Read and base64-encode an image file."""
        path = Path(image_path)
        if not path.exists():
            raise FileNotFoundError(f"Image not found: {image_path}")
        with open(path, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")

    def _get_mime_type(self, image_path: str) -> str:
        """Determine MIME type from file extension."""
        ext = Path(image_path).suffix.lower()
        mime_map = {
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".dcm": "application/dicom",
            ".nii": "application/octet-stream",
            ".nii.gz": "application/octet-stream",
        }
        return mime_map.get(ext, "image/png")

    def _build_vlm_payload(
        self,
        image_path: str,
        prompt: str,
        temperature: float = 0.1,
        max_tokens: int = 1024,
    ) -> Dict[str, Any]:
        """Build an OpenAI-compatible chat completion payload with image."""
        image_b64 = self._encode_image(image_path)
        mime_type = self._get_mime_type(image_path)

        return {
            "model": "vila-m3",
            "messages": [
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "image_url",
                            "image_url": {
                                "url": f"data:{mime_type};base64,{image_b64}"
                            },
                        },
                        {
                            "type": "text",
                            "text": prompt,
                        },
                    ],
                }
            ],
            "temperature": temperature,
            "max_tokens": max_tokens,
        }

    def _build_cloud_vlm_payload(
        self,
        image_path: str,
        prompt: str,
        temperature: float = 0.1,
        max_tokens: int = 1024,
    ) -> Dict[str, Any]:
        """Build payload for NVIDIA cloud Llama 3.2 Vision endpoint.

        Uses OpenAI vision API compatible format with base64 image.
        """
        image_b64 = self._encode_image(image_path)
        mime_type = self._get_mime_type(image_path)

        return {
            "model": self.cloud_vlm_model,
            "messages": [
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "text",
                            "text": prompt,
                        },
                        {
                            "type": "image_url",
                            "image_url": {
                                "url": f"data:{mime_type};base64,{image_b64}"
                            },
                        },
                    ],
                }
            ],
            "max_tokens": max_tokens,
            "temperature": temperature,
        }

    def _cloud_analyze(
        self,
        image_path: str,
        prompt: str,
        temperature: float = 0.1,
        max_tokens: int = 1024,
    ) -> Optional[Dict[str, Any]]:
        """Call NVIDIA cloud VLM endpoint for image analysis.

        Args:
            image_path: Path to the image file.
            prompt: Text prompt for the VLM.
            temperature: Sampling temperature.
            max_tokens: Maximum response tokens.

        Returns:
            OpenAI-compatible response dict, or None on failure.
        """
        if not self.nvidia_api_key:
            return None

        try:
            payload = self._build_cloud_vlm_payload(
                image_path, prompt, temperature, max_tokens
            )

            resp = requests.post(
                f"{self.cloud_url}/chat/completions",
                json=payload,
                headers={
                    "Authorization": f"Bearer {self.nvidia_api_key}",
                    "Content-Type": "application/json",
                },
                timeout=120,
            )
            resp.raise_for_status()
            result = resp.json()

            content = ""
            if "choices" in result and len(result["choices"]) > 0:
                content = result["choices"][0].get("message", {}).get("content", "")

            logger.info(
                f"Cloud NIM VLM ({self.cloud_vlm_model}) generated {len(content)} chars"
            )
            return result
        except Exception as e:
            logger.error(f"Cloud NIM VLM request failed: {e}")
            return None

    def _invoke_or_mock(
        self, endpoint: str, payload: Dict, timeout: int = 120, **mock_kwargs
    ) -> Any:
        """Try local NIM, then cloud NIM, then fall back to mock.

        Extends the base class to add cloud NIM as a middle fallback.
        """
        # Attempt 1: Local NIM
        if self.is_available():
            try:
                return self._request(endpoint, payload, timeout)
            except Exception as e:
                logger.error(f"NIM {self.service_name} request failed: {e}")

        # Attempt 2: NVIDIA Cloud NIM
        if self.nvidia_api_key:
            image_path = mock_kwargs.get("image_path")
            question = mock_kwargs.get("question", "")
            mode = mock_kwargs.get("mode", "vqa")

            if image_path and Path(image_path).exists():
                # Build the prompt based on mode
                if mode == "report":
                    prompt = (
                        "You are an expert radiologist. Generate a structured "
                        "radiology report for this medical image. Include sections: "
                        "TECHNIQUE, COMPARISON, FINDINGS, IMPRESSION."
                    )
                    findings_context = mock_kwargs.get("findings_context", "")
                    if findings_context:
                        prompt += f"\n\nAdditional clinical context:\n{findings_context}"
                elif mode == "classify":
                    labels = mock_kwargs.get("labels", [])
                    labels_str = ", ".join(f'"{lbl}"' for lbl in labels)
                    prompt = (
                        f"Classify this medical image into one of: [{labels_str}]. "
                        f"Provide confidence scores 0.0-1.0 for each in JSON format."
                    )
                else:
                    prompt = question if question else "Describe this medical image."

                cloud_result = self._cloud_analyze(image_path, prompt)
                if cloud_result is not None:
                    return cloud_result

        # Attempt 3: Mock
        if self.mock_enabled:
            logger.info(f"Using mock for {self.service_name} (local and cloud unavailable)")
            return self._mock_response(**mock_kwargs)

        raise ConnectionError(
            f"NIM {self.service_name} unavailable: local NIM unreachable, "
            f"cloud NIM failed, and mock disabled"
        )

    def analyze_image(
        self,
        image_path: str,
        question: str,
    ) -> VLMResponse:
        """Visual question answering on a medical image.

        Args:
            image_path: Path to the medical image (PNG, JPEG, DICOM, NIfTI).
            question: Clinical question about the image.

        Returns:
            VLMResponse with answer text, confidence, and metadata.
        """
        logger.info(f"VQA on {image_path}: {question[:80]}...")

        # Build payload only when not using mock (avoids file-not-found in mock mode)
        if self.is_available():
            payload = self._build_vlm_payload(image_path, question)
        else:
            payload = {}

        result = self._invoke_or_mock(
            endpoint="/v1/chat/completions",
            payload=payload,
            timeout=120,
            image_path=image_path,
            question=question,
            mode="vqa",
        )

        # Parse OpenAI-compatible response into VLMResponse
        if isinstance(result, dict) and not isinstance(result, VLMResponse):
            answer = ""
            if "choices" in result and len(result["choices"]) > 0:
                answer = result["choices"][0].get("message", {}).get("content", "")

            # Determine the model name for tracking
            model_name = result.get("model", "")
            if not model_name:
                model_name = "vila-m3"
            is_cloud = self.cloud_vlm_model in str(model_name)

            findings = [s.strip() + "." for s in answer.split(".") if s.strip()]
            return VLMResponse(
                answer=answer,
                findings=findings[:5],
                confidence=result.get("confidence", 0.0),
                model=f"cloud:{self.cloud_vlm_model}" if is_cloud else "vila-m3",
                is_mock=False,
            )

        return result

    def generate_report(
        self,
        image_path: str,
        findings_context: str = "",
    ) -> str:
        """Generate a structured radiology report from a medical image.

        Args:
            image_path: Path to the medical image.
            findings_context: Optional prior findings or clinical context
                to incorporate into the report.

        Returns:
            Structured radiology report as a string.
        """
        logger.info(f"Generating radiology report for: {image_path}")

        prompt = (
            "You are an expert radiologist. Generate a structured radiology report "
            "for this medical image. Include the following sections:\n"
            "1. TECHNIQUE\n"
            "2. COMPARISON\n"
            "3. FINDINGS\n"
            "4. IMPRESSION\n\n"
            "Be specific about anatomical structures, measurements, and any "
            "abnormalities observed."
        )
        if findings_context:
            prompt += f"\n\nAdditional clinical context:\n{findings_context}"

        if self.is_available():
            payload = self._build_vlm_payload(
                image_path, prompt, temperature=0.1, max_tokens=2048
            )
        else:
            payload = {}

        result = self._invoke_or_mock(
            endpoint="/v1/chat/completions",
            payload=payload,
            timeout=180,
            image_path=image_path,
            findings_context=findings_context,
            mode="report",
        )

        # Extract report text
        if isinstance(result, dict):
            if "choices" in result and len(result["choices"]) > 0:
                return result["choices"][0].get("message", {}).get("content", "")
            return result.get("report", "")
        if isinstance(result, VLMResponse):
            return result.answer

        return str(result)

    def classify(
        self,
        image_path: str,
        labels: List[str],
    ) -> Dict[str, float]:
        """Zero-shot classification of a medical image.

        Args:
            image_path: Path to the medical image.
            labels: List of candidate classification labels.

        Returns:
            Dictionary mapping each label to its probability score.
        """
        logger.info(
            f"Zero-shot classification with {len(labels)} labels: {image_path}"
        )

        labels_str = ", ".join(f'"{lbl}"' for lbl in labels)
        prompt = (
            f"Classify this medical image into exactly one of the following "
            f"categories: [{labels_str}]. For each category, provide a "
            f"confidence score between 0.0 and 1.0. Respond in JSON format "
            f'like: {{"label1": 0.85, "label2": 0.10, ...}}'
        )

        if self.is_available():
            payload = self._build_vlm_payload(
                image_path, prompt, temperature=0.0, max_tokens=512
            )
        else:
            payload = {}

        result = self._invoke_or_mock(
            endpoint="/v1/chat/completions",
            payload=payload,
            timeout=120,
            image_path=image_path,
            labels=labels,
            mode="classify",
        )

        # Parse classification scores from response
        content = ""
        if isinstance(result, VLMResponse):
            content = result.answer
        elif isinstance(result, dict):
            if "choices" in result and len(result["choices"]) > 0:
                content = result["choices"][0].get("message", {}).get("content", "")
            elif all(lbl in result for lbl in labels):
                return {lbl: float(result[lbl]) for lbl in labels}

        if content:
            try:
                import json
                scores = json.loads(content)
                if isinstance(scores, dict):
                    return {k: float(v) for k, v in scores.items()}
            except (json.JSONDecodeError, ValueError):
                logger.warning(
                    f"Could not parse classification JSON from VLM: {content[:200]}"
                )

        # Fallback: return uniform scores
        logger.warning("Could not extract classification scores, returning uniform")
        uniform = 1.0 / len(labels) if labels else 0.0
        return {lbl: uniform for lbl in labels}

    def get_status(self) -> str:
        """Return service status string.

        Returns "available" for local NIM, "cloud" for NVIDIA cloud NIM,
        "mock" for mock mode, or "unavailable".
        """
        if self.is_available():
            return "available"
        if self.nvidia_api_key and self.is_cloud_available():
            return "cloud"
        if self.mock_enabled:
            return "mock"
        return "unavailable"

    def _mock_response(self, **kwargs) -> VLMResponse:
        """Return a realistic mock VLMResponse.

        Generates plausible radiology findings based on the mode
        (VQA, report generation, or classification).
        """
        mode = kwargs.get("mode", "vqa")
        image_path = kwargs.get("image_path", "unknown.png")
        question = kwargs.get("question", "")

        if mode == "report":
            answer = (
                "TECHNIQUE: CT of the chest without contrast.\n\n"
                "COMPARISON: No prior studies available.\n\n"
                "FINDINGS:\n"
                "- Lungs: Clear bilaterally. No focal consolidation, "
                "pleural effusion, or pneumothorax.\n"
                "- Heart: Normal size. No pericardial effusion.\n"
                "- Mediastinum: No lymphadenopathy. Normal vascular structures.\n"
                "- Bones: No acute fracture or destructive lesion.\n"
                "- Upper abdomen: Limited evaluation; visualized portions "
                "appear unremarkable.\n\n"
                "IMPRESSION:\n"
                "1. No acute cardiopulmonary abnormality.\n"
                "2. Normal chest CT examination."
            )
        elif mode == "classify":
            labels = kwargs.get("labels", ["normal", "abnormal"])
            # Generate mock classification scores
            import json
            primary_idx = random.randint(0, len(labels) - 1)
            scores = {}
            remaining = 1.0
            for i, lbl in enumerate(labels):
                if i == primary_idx:
                    scores[lbl] = round(random.uniform(0.65, 0.90), 3)
                    remaining -= scores[lbl]
                else:
                    score = round(remaining / (len(labels) - i), 3)
                    scores[lbl] = max(0.01, score)
                    remaining -= scores[lbl]

            answer = json.dumps(scores, indent=2)
        else:
            # VQA mode
            findings = [
                "The image shows normal lung parenchyma with no evidence "
                "of consolidation or ground-glass opacities.",
                "The cardiac silhouette is within normal limits.",
                "No pleural effusion or pneumothorax is identified.",
                "The mediastinal structures appear normal.",
                "The visualized osseous structures are intact without "
                "acute fracture.",
            ]
            selected = random.sample(findings, min(3, len(findings)))
            answer = " ".join(selected)

        logger.info(f"Mock VILA-M3 response generated ({mode} mode)")

        # Build findings list from the answer text
        findings_list = [s.strip() for s in answer.split("\n") if s.strip() and s.strip() != "-"]
        if not findings_list:
            findings_list = [s.strip() + "." for s in answer.split(".") if s.strip()]

        return VLMResponse(
            answer=answer,
            findings=findings_list[:5],
            confidence=round(random.uniform(0.70, 0.95), 3),
            model="vila-m3-mock",
            is_mock=True,
        )
