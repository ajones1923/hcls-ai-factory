"""Live DICOM analysis API routes for the Clinical Imaging Engine.

Provides endpoints to upload DICOM files, run live AI inference,
and analyze pre-loaded sample images. All results are real inference
(is_mock=False), distinct from the canned demo workflows.

Author: Adam Jones
Date: April 2026
"""

import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional

import torch
from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from loguru import logger

from src.dicom_analyzer import DICOMAnalyzer

# =====================================================================
# Router
# =====================================================================

router = APIRouter(prefix="/analyze", tags=["Live Analysis"])

# Singleton analyzer (lazy init with GPU if available)
_analyzer: Optional[DICOMAnalyzer] = None

# Project root for sample data paths
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def _get_analyzer() -> DICOMAnalyzer:
    """Get or create the singleton DICOMAnalyzer."""
    global _analyzer
    if _analyzer is None:
        models_dir = str(_PROJECT_ROOT / "data" / "models")
        _analyzer = DICOMAnalyzer(models_dir=models_dir, use_gpu=True)
    return _analyzer


# Pre-loaded sample DICOM paths
_SAMPLE_PATHS = {
    "sample_cxr": "data/sample_images/sample_cxr.dcm",
    "cxr_maria_santos": "data/demo/dicom/cxr_maria_santos/cxr_001.dcm",
    "ct_chest_maria_santos": "data/demo/dicom/ct_chest_maria_santos",
}


# =====================================================================
# Endpoints
# =====================================================================


@router.post("/upload")
async def analyze_uploaded_dicom(
    file: UploadFile = File(...),
    workflow: Optional[str] = Form(None),
):
    """Upload a single DICOM file for live AI analysis.

    If workflow is not specified, the analyzer auto-detects the appropriate
    workflow based on DICOM modality and body part headers.

    Returns real AI analysis results (is_mock=False).
    """
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file uploaded")

    # Save to temp location
    suffix = ".dcm" if not file.filename.lower().endswith(".dcm") else ""
    tmp_dir = tempfile.mkdtemp(prefix="dicom_upload_")
    tmp_path = os.path.join(tmp_dir, file.filename + suffix)

    try:
        content = await file.read()
        with open(tmp_path, "wb") as f:
            f.write(content)
        logger.info(f"Uploaded DICOM saved: {tmp_path} ({len(content)} bytes)")

        analyzer = _get_analyzer()
        result = analyzer.analyze_dicom(tmp_path, workflow_name=workflow)
        return result

    except Exception as e:
        logger.error(f"Upload analysis failed: {e}")
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")
    finally:
        # Clean up temp files
        try:
            shutil.rmtree(tmp_dir, ignore_errors=True)
        except Exception:
            pass


@router.post("/upload-series")
async def analyze_uploaded_series(
    files: list[UploadFile] = File(...),
    workflow: Optional[str] = Form(None),
):
    """Upload multiple DICOM slices (CT/MRI series) for live analysis.

    All files are saved to a temporary directory and analyzed as a volume.
    """
    if not files:
        raise HTTPException(status_code=400, detail="No files uploaded")

    tmp_dir = tempfile.mkdtemp(prefix="dicom_series_")

    try:
        # Save all slices
        for i, f in enumerate(files):
            fname = f.filename or f"slice_{i:04d}.dcm"
            fpath = os.path.join(tmp_dir, fname)
            content = await f.read()
            with open(fpath, "wb") as out:
                out.write(content)

        logger.info(f"Uploaded {len(files)} DICOM slices to {tmp_dir}")

        analyzer = _get_analyzer()
        result = analyzer.analyze_dicom(tmp_dir, workflow_name=workflow)
        return result

    except Exception as e:
        logger.error(f"Series analysis failed: {e}")
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")
    finally:
        try:
            shutil.rmtree(tmp_dir, ignore_errors=True)
        except Exception:
            pass


@router.post("/analyze-sample/{sample_name}")
async def analyze_sample(sample_name: str):
    """Run live analysis on a pre-loaded sample DICOM.

    Available samples:
    - sample_cxr: Sample chest X-ray DICOM
    - cxr_maria_santos: Demo CXR for Maria Santos case
    - ct_chest_maria_santos: 64-slice CT chest series (Maria Santos)
    """
    if sample_name not in _SAMPLE_PATHS:
        raise HTTPException(
            status_code=404,
            detail=f"Sample '{sample_name}' not found. "
            f"Available: {list(_SAMPLE_PATHS.keys())}",
        )

    rel_path = _SAMPLE_PATHS[sample_name]
    abs_path = str(_PROJECT_ROOT / rel_path)

    if not os.path.exists(abs_path):
        raise HTTPException(
            status_code=404,
            detail=f"Sample file not found on disk: {rel_path}",
        )

    logger.info(f"Analyzing sample '{sample_name}': {abs_path}")

    analyzer = _get_analyzer()

    try:
        result = analyzer.analyze_dicom(abs_path)
        result["sample_name"] = sample_name
        return result
    except Exception as e:
        logger.error(f"Sample analysis failed: {e}")
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")


@router.get("/status")
async def analysis_status():
    """Check if live analysis is available (GPU, models loaded, etc.)."""
    analyzer = _get_analyzer()
    status = analyzer.get_status()

    # Add available samples
    available_samples = {}
    for name, rel_path in _SAMPLE_PATHS.items():
        abs_path = str(_PROJECT_ROOT / rel_path)
        available_samples[name] = {
            "path": rel_path,
            "exists": os.path.exists(abs_path),
            "is_series": os.path.isdir(abs_path),
        }
    status["available_samples"] = available_samples

    return status


@router.get("/samples")
async def list_samples():
    """List all available pre-loaded sample DICOMs with metadata."""
    samples = []
    for name, rel_path in _SAMPLE_PATHS.items():
        abs_path = str(_PROJECT_ROOT / rel_path)
        exists = os.path.exists(abs_path)
        is_series = os.path.isdir(abs_path) if exists else False

        sample_info = {
            "name": name,
            "path": rel_path,
            "exists": exists,
            "is_series": is_series,
        }

        if exists and is_series:
            dcm_count = len([
                f for f in os.listdir(abs_path)
                if f.lower().endswith(".dcm") or "." not in f
            ])
            sample_info["slice_count"] = dcm_count

        samples.append(sample_info)

    return {"samples": samples}
