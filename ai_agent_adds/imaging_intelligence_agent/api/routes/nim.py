"""NIM service endpoints for Imaging Intelligence Agent.

Provides health checks and proxy endpoints for NVIDIA NIM services:
  - VISTA-3D: 3D medical image segmentation
  - MAISI: Synthetic CT generation via latent diffusion
  - VILA-M3: Vision-language model for radiology

Author: Adam Jones
Date: February 2026
"""

import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, UploadFile, File
from loguru import logger
from pydantic import BaseModel, Field

# =====================================================================
# Request / Response Models
# =====================================================================


class NIMServiceStatus(BaseModel):
    """Status of a single NIM service."""
    name: str
    status: str = Field(..., description="available, cloud, anthropic, mock, or unavailable")
    url: str = ""


class NIMStatusResponse(BaseModel):
    """Aggregated NIM service health status."""
    services: List[NIMServiceStatus]
    available_count: int
    mock_count: int
    unavailable_count: int


class SegmentRequest(BaseModel):
    """VISTA-3D segmentation request."""
    input_path: str = Field("", description="Path to NIfTI volume (empty for mock)")
    target_classes: Optional[List[str]] = Field(
        None, description="Specific anatomy classes to segment"
    )


class SegmentResponse(BaseModel):
    """VISTA-3D segmentation response."""
    classes_detected: List[str]
    volumes: Dict[str, float] = Field(default_factory=dict, description="Class -> volume in mL")
    num_classes: int
    inference_time_ms: float
    model: str = "vista3d"
    is_mock: bool


class GenerateRequest(BaseModel):
    """MAISI synthetic CT generation request."""
    body_region: str = Field("chest", description="Target body region")
    resolution: str = Field("512x512x512", description="Output volume resolution")
    num_classes: int = Field(104, ge=1, le=200, description="Number of annotation classes")


class GenerateResponse(BaseModel):
    """MAISI synthetic CT generation response."""
    resolution: str
    body_region: str
    num_classes_annotated: int
    generation_time_ms: float
    model: str = "maisi"
    is_mock: bool


class AnalyzeRequest(BaseModel):
    """VILA-M3 image analysis request."""
    question: str = Field(
        "Describe the findings in this image.",
        description="Question to ask about the image",
    )
    input_path: str = Field("", description="Path to image (empty for mock)")


class AnalyzeResponse(BaseModel):
    """VILA-M3 image analysis response."""
    answer: str
    findings: List[str]
    confidence: float
    inference_time_ms: float
    model: str = "vila_m3"
    is_mock: bool


# =====================================================================
# Router
# =====================================================================

router = APIRouter()


def _get_state():
    """Import and access the shared application state from api.main."""
    from api.main import _state
    return _state


def _get_nim_manager():
    """Get NIM service manager or raise 503."""
    state = _get_state()
    nim_manager = state.get("nim_manager")
    if nim_manager is None:
        raise HTTPException(status_code=503, detail="NIM service manager not initialized")
    return nim_manager


# =====================================================================
# Endpoints
# =====================================================================


@router.get("/status", response_model=NIMStatusResponse)
async def nim_status():
    """Health check all NIM services.

    Returns the status of each NIM service (VISTA-3D, MAISI, VILA-M3, LLM)
    along with aggregate counts.
    """
    nim_manager = _get_nim_manager()
    settings = _get_state().get("settings")

    try:
        status_dict = nim_manager.check_all_services()
    except Exception as e:
        logger.error(f"NIM status check failed: {e}")
        raise HTTPException(status_code=500, detail=f"NIM status check failed: {e}")

    cloud_url = getattr(settings, "NIM_CLOUD_URL", "")
    nim_urls = {
        "vista3d": getattr(settings, "NIM_VISTA3D_URL", ""),
        "maisi": getattr(settings, "NIM_MAISI_URL", ""),
        "vila_m3": getattr(settings, "NIM_VILAM3_URL", ""),
        "llm": getattr(settings, "NIM_LLM_URL", ""),
    }

    services = []
    for name, status in status_dict.items():
        # Show cloud URL when the service is using cloud mode
        url = cloud_url if status == "cloud" else nim_urls.get(name, "")
        services.append(NIMServiceStatus(
            name=name,
            status=status,
            url=url,
        ))

    available = sum(1 for s in services if s.status in ("available", "cloud", "anthropic"))
    mock = sum(1 for s in services if s.status == "mock")
    unavailable = sum(1 for s in services if s.status == "unavailable")

    return NIMStatusResponse(
        services=services,
        available_count=available,
        mock_count=mock,
        unavailable_count=unavailable,
    )


@router.post("/vista3d/segment", response_model=SegmentResponse)
async def vista3d_segment(request: SegmentRequest):
    """Run VISTA-3D 3D segmentation on a medical volume.

    Supports 132 anatomical classes with zero-shot and interactive prompting.
    Falls back to mock if the NIM endpoint is unavailable.
    """
    nim_manager = _get_nim_manager()
    start = time.time()

    try:
        client = nim_manager.vista3d
        result = client.segment(
            ct_volume_path=request.input_path or "",
            classes=request.target_classes if request.target_classes else None,
        )

        return SegmentResponse(
            classes_detected=result.classes_detected,
            volumes=result.volumes,
            num_classes=result.num_classes,
            inference_time_ms=(time.time() - start) * 1000,
            model=result.model,
            is_mock=result.is_mock,
        )
    except Exception as e:
        logger.error(f"VISTA-3D segmentation failed: {e}")
        raise HTTPException(status_code=500, detail=f"Segmentation failed: {e}")


@router.post("/maisi/generate", response_model=GenerateResponse)
async def maisi_generate(request: GenerateRequest):
    """Generate synthetic CT volume with MAISI latent diffusion model.

    Produces high-resolution synthetic CT data with paired segmentation masks
    for up to 127 anatomical classes.
    Falls back to mock if the NIM endpoint is unavailable.
    """
    nim_manager = _get_nim_manager()
    start = time.time()

    try:
        client = nim_manager.maisi
        result = client.generate(
            body_region=request.body_region,
            resolution=request.resolution,
        )

        return GenerateResponse(
            resolution=result.resolution,
            body_region=result.body_region,
            num_classes_annotated=result.num_classes_annotated,
            generation_time_ms=(time.time() - start) * 1000,
            model=result.model,
            is_mock=result.is_mock,
        )
    except Exception as e:
        logger.error(f"MAISI generation failed: {e}")
        raise HTTPException(status_code=500, detail=f"Synthetic CT generation failed: {e}")


@router.post("/vilam3/analyze", response_model=AnalyzeResponse)
async def vilam3_analyze(request: AnalyzeRequest):
    """Analyze a medical image using VILA-M3 vision-language model.

    Provides natural language understanding of radiology images,
    combining visual features with medical domain knowledge.
    Falls back to mock if the NIM endpoint is unavailable.
    """
    nim_manager = _get_nim_manager()
    start = time.time()

    try:
        client = nim_manager.vilam3
        result = client.analyze_image(
            image_path=request.input_path or "",
            question=request.question,
        )

        findings = result.findings if result.findings else [
            s.strip() for s in result.answer.split(".") if s.strip()
        ]
        return AnalyzeResponse(
            answer=result.answer,
            findings=findings[:5],
            confidence=result.confidence,
            inference_time_ms=(time.time() - start) * 1000,
            model=result.model,
            is_mock=result.is_mock,
        )
    except Exception as e:
        logger.error(f"VILA-M3 analysis failed: {e}")
        raise HTTPException(status_code=500, detail=f"Image analysis failed: {e}")
