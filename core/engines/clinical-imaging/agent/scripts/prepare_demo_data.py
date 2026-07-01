#!/usr/bin/env python3
"""Prepare all demo data for Clinical Imaging Engine demonstrations.

Generates:
1. Multi-slice DICOM CT chest study (Maria Santos case) from sample NIfTI
2. Pre-computed workflow results for all 9 workflows (JSON)
3. Pre-computed radiomics features (JSON)
4. Sample radiology reports (text files)
5. Pre-computed report NLP results (JSON)
6. Pre-computed cross-modal genomic trigger results (JSON)
7. Pre-computed analytics demo data (JSON)
8. Pre-computed DICOM SR files
9. Demo target hypothesis for Engine 3 handoff (JSON)
10. CXR DICOM study
11. Dose tracking records

Usage:
    python scripts/prepare_demo_data.py
    python scripts/prepare_demo_data.py --output-dir data/demo
    python scripts/prepare_demo_data.py --patient-name "Maria Santos" --patient-id "MS2026"

Author: Adam Jones
Date: April 2026
"""

import argparse
import json
import os
import random
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

SAMPLE_DIR = PROJECT_ROOT / "data" / "sample_images"

DEFAULT_PATIENT_NAME = "Maria Santos"
DEFAULT_PATIENT_ID = "MS2026"
DEFAULT_NUM_CT_SLICES = 64

# Maria Santos demographics
MARIA_BIRTH_DATE = "19680315"
MARIA_SEX = "F"
MARIA_WEIGHT = "65"
MARIA_AGE = 58

TODAY = datetime.now()
STUDY_DATE = TODAY.strftime("%Y%m%d")
STUDY_TIME = TODAY.strftime("%H%M%S")


# ═══════════════════════════════════════════════════════════════════════
# 1. DICOM CT CHEST STUDY
# ═══════════════════════════════════════════════════════════════════════


def generate_ct_chest_dicom_study(
    output_dir: Path,
    patient_name: str = DEFAULT_PATIENT_NAME,
    patient_id: str = DEFAULT_PATIENT_ID,
    num_slices: int = DEFAULT_NUM_CT_SLICES,
) -> Path:
    """Generate a synthetic multi-slice CT chest DICOM study.

    Creates proper DICOM files with:
    - Correct UIDs (StudyInstanceUID, SeriesInstanceUID, SOPInstanceUID per slice)
    - Patient demographics: Maria Santos, 58F, 65kg
    - Study: CT Chest with Contrast, today's date
    - Series: Axial, 1.25mm slice thickness, 512x512 matrix
    - Institution: Community Hospital, Rural NM
    - Pixel data: From sample_ct_chest.nii.gz if available, else synthetic noise

    Args:
        output_dir: Root output directory (will create dicom/ct_chest_maria_santos/).
        patient_name: Patient name for DICOM tags.
        patient_id: Patient ID for DICOM tags.
        num_slices: Number of CT slices to generate.

    Returns:
        Path to the directory containing generated DICOM files.
    """
    import pydicom
    from pydicom.uid import generate_uid, ExplicitVRLittleEndian

    dicom_dir = output_dir / "dicom" / "ct_chest_maria_santos"
    dicom_dir.mkdir(parents=True, exist_ok=True)

    # Format patient name as DICOM PN (Last^First)
    name_parts = patient_name.split()
    dicom_patient_name = f"{name_parts[-1].upper()}^{name_parts[0].upper()}" if len(name_parts) >= 2 else patient_name.upper()

    # Try loading NIfTI volume
    nifti_path = SAMPLE_DIR / "sample_ct_chest.nii.gz"
    pixel_source = "synthetic"
    rows, cols = 256, 256

    if nifti_path.exists():
        try:
            import nibabel as nib

            img = nib.load(str(nifti_path))
            vol = img.get_fdata()
            rows, cols = vol.shape[0], vol.shape[1]
            total_z = vol.shape[2]

            # Select evenly-spaced slices
            step = max(1, total_z // num_slices)
            slice_indices = list(range(0, total_z, step))[:num_slices]
            num_slices = len(slice_indices)

            pixel_source = "nifti"
            logger.info(f"Loaded NIfTI volume: {vol.shape}, extracting {num_slices} slices")
        except Exception as e:
            logger.warning(f"Failed to load NIfTI ({e}), using synthetic data")
            pixel_source = "synthetic"
    else:
        logger.info("NIfTI sample not found, generating synthetic CT data")

    # Generate shared UIDs
    study_uid = generate_uid()
    series_uid = generate_uid()
    frame_of_ref_uid = generate_uid()

    logger.info(f"Generating {num_slices} DICOM CT slices -> {dicom_dir}")

    for i in range(num_slices):
        sop_uid = generate_uid()

        # File meta
        file_meta = pydicom.dataset.FileMetaDataset()
        file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"  # CT Image Storage
        file_meta.MediaStorageSOPInstanceUID = sop_uid

        ds = pydicom.dataset.FileDataset(
            filename_or_obj="",
            dataset={},
            file_meta=file_meta,
            preamble=b"\x00" * 128,
        )

        # Patient module
        ds.PatientName = dicom_patient_name
        ds.PatientID = patient_id
        ds.PatientBirthDate = MARIA_BIRTH_DATE
        ds.PatientSex = MARIA_SEX
        ds.PatientWeight = MARIA_WEIGHT
        ds.PatientAge = f"{MARIA_AGE:03d}Y"

        # Study module
        ds.StudyInstanceUID = study_uid
        ds.StudyDate = STUDY_DATE
        ds.StudyTime = STUDY_TIME
        ds.StudyDescription = "CT CHEST WITH CONTRAST"
        ds.StudyID = "1"
        ds.AccessionNumber = f"ACC{patient_id}"
        ds.ReferringPhysicianName = "RODRIGUEZ^ELENA^MD"
        ds.InstitutionName = "Community Hospital - Rural NM"
        ds.Manufacturer = "HCLS AI Factory Demo"
        ds.ManufacturerModelName = "Clinical Imaging Engine v2.0"

        # Series module
        ds.SeriesInstanceUID = series_uid
        ds.SeriesNumber = 3
        ds.SeriesDescription = "AXIAL 1.25MM"
        ds.Modality = "CT"
        ds.FrameOfReferenceUID = frame_of_ref_uid
        ds.BodyPartExamined = "CHEST"
        ds.ProtocolName = "CT CHEST W CONTRAST"

        # Instance
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        ds.SOPInstanceUID = sop_uid
        ds.InstanceNumber = i + 1
        slice_location = float(i * 1.25)
        ds.ImagePositionPatient = [0.0, 0.0, slice_location]
        ds.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        ds.SliceLocation = slice_location
        ds.SliceThickness = "1.25"
        ds.PixelSpacing = [0.703125, 0.703125]

        # CT-specific
        ds.KVP = "120"
        ds.XRayTubeCurrent = 200
        ds.ConvolutionKernel = "STANDARD"
        ds.ReconstructionDiameter = 360.0

        # Pixel data
        if pixel_source == "nifti":
            z_idx = slice_indices[i]
            raw_slice = vol[:, :, z_idx].astype(np.float64)
            # Map to signed 16-bit with HU offset
            # Apply RescaleIntercept=-1024, RescaleSlope=1
            # So stored = (HU + 1024)
            pixel_data = np.clip(raw_slice + 1024, 0, 4095).astype(np.int16)
        else:
            # Synthetic CT-like data
            pixel_data = _generate_synthetic_ct_slice(rows, cols, i, num_slices)

        ds.Rows = pixel_data.shape[0]
        ds.Columns = pixel_data.shape[1]
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1  # signed
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.RescaleIntercept = "-1024"
        ds.RescaleSlope = "1"
        ds.WindowCenter = "40"
        ds.WindowWidth = "400"
        ds.PixelData = pixel_data.astype(np.int16).tobytes()

        filepath = dicom_dir / f"slice_{i + 1:03d}.dcm"
        ds.save_as(str(filepath), write_like_original=False)

    logger.info(f"Generated {num_slices} DICOM CT slices in {dicom_dir}")
    return dicom_dir


def _generate_synthetic_ct_slice(
    rows: int, cols: int, slice_idx: int, total_slices: int
) -> np.ndarray:
    """Generate a synthetic CT-like slice with a circular body, lungs, and optional nodule.

    Values are stored pixel values (HU + 1024):
    - Air (outside body): 0 (HU -1024)
    - Lung parenchyma: ~124 (HU ~-900)
    - Soft tissue: ~1064 (HU ~40)
    - Bone (spine): ~1324 (HU ~300)
    - Nodule: ~1084 (HU ~60)

    Args:
        rows: Image height.
        cols: Image width.
        slice_idx: Current slice index.
        total_slices: Total number of slices.

    Returns:
        np.ndarray of shape (rows, cols), dtype int16.
    """
    rng = np.random.RandomState(42 + slice_idx)
    img = np.zeros((rows, cols), dtype=np.float64)

    cy, cx = rows // 2, cols // 2
    yy, xx = np.ogrid[:rows, :cols]

    # Body ellipse
    body_ry, body_rx = int(rows * 0.38), int(cols * 0.42)
    body_mask = ((yy - cy) / body_ry) ** 2 + ((xx - cx) / body_rx) ** 2 <= 1.0
    img[body_mask] = 1064  # soft tissue ~ 40 HU

    # Left lung
    lung_l_cy, lung_l_cx = cy - int(rows * 0.02), cx - int(cols * 0.15)
    lung_l_ry, lung_l_rx = int(rows * 0.25), int(cols * 0.12)
    lung_l_mask = ((yy - lung_l_cy) / lung_l_ry) ** 2 + ((xx - lung_l_cx) / lung_l_rx) ** 2 <= 1.0
    img[lung_l_mask] = 124  # lung ~ -900 HU

    # Right lung
    lung_r_cy, lung_r_cx = cy - int(rows * 0.02), cx + int(cols * 0.15)
    lung_r_ry, lung_r_rx = int(rows * 0.25), int(cols * 0.12)
    lung_r_mask = ((yy - lung_r_cy) / lung_r_ry) ** 2 + ((xx - lung_r_cx) / lung_r_rx) ** 2 <= 1.0
    img[lung_r_mask] = 124

    # Spine (posterior midline)
    spine_cy = cy + int(rows * 0.28)
    spine_r = int(rows * 0.06)
    spine_mask = (yy - spine_cy) ** 2 + (xx - cx) ** 2 <= spine_r ** 2
    img[spine_mask] = 1324  # bone ~ 300 HU

    # Nodule in right upper lobe region (slices 15-25 of a 64-slice study)
    frac = slice_idx / max(total_slices - 1, 1)
    if 0.2 <= frac <= 0.4:
        nodule_cy = cy - int(rows * 0.12)
        nodule_cx = cx + int(cols * 0.10)
        nodule_r = max(2, int(rows * 0.015))
        nodule_mask = (yy - nodule_cy) ** 2 + (xx - nodule_cx) ** 2 <= nodule_r ** 2
        img[nodule_mask] = 1084  # nodule ~ 60 HU

    # Add realistic noise
    noise = rng.normal(0, 8, (rows, cols))
    img += noise
    img = np.clip(img, 0, 4095)

    return img.astype(np.int16)


# ═══════════════════════════════════════════════════════════════════════
# 2. CXR DICOM STUDY
# ═══════════════════════════════════════════════════════════════════════


def generate_cxr_dicom(
    output_dir: Path,
    patient_name: str = DEFAULT_PATIENT_NAME,
    patient_id: str = DEFAULT_PATIENT_ID,
) -> Path:
    """Generate a CXR DICOM from sample_cxr.dcm or sample PNG.

    Args:
        output_dir: Root output directory.
        patient_name: Patient name.
        patient_id: Patient ID.

    Returns:
        Path to the generated CXR DICOM file.
    """
    import pydicom
    from pydicom.uid import generate_uid, ExplicitVRLittleEndian

    cxr_dir = output_dir / "dicom" / "cxr_maria_santos"
    cxr_dir.mkdir(parents=True, exist_ok=True)

    name_parts = patient_name.split()
    dicom_patient_name = f"{name_parts[-1].upper()}^{name_parts[0].upper()}" if len(name_parts) >= 2 else patient_name.upper()

    # Try copying from existing sample DICOM, re-tagging patient info
    sample_cxr = SAMPLE_DIR / "sample_cxr.dcm"

    if sample_cxr.exists():
        try:
            ds = pydicom.dcmread(str(sample_cxr))
            ds.PatientName = dicom_patient_name
            ds.PatientID = patient_id
            ds.PatientBirthDate = MARIA_BIRTH_DATE
            ds.PatientSex = MARIA_SEX
            ds.PatientWeight = MARIA_WEIGHT
            ds.PatientAge = f"{MARIA_AGE:03d}Y"
            ds.StudyDescription = "CHEST PA AND LATERAL"
            ds.SeriesDescription = "PA UPRIGHT"
            ds.InstitutionName = "Community Hospital - Rural NM"
            ds.Manufacturer = "HCLS AI Factory Demo"
            ds.ReferringPhysicianName = "RODRIGUEZ^ELENA^MD"
            ds.StudyDate = STUDY_DATE
            ds.StudyTime = STUDY_TIME

            out_path = cxr_dir / "cxr_001.dcm"
            ds.save_as(str(out_path), write_like_original=False)
            logger.info(f"Generated CXR DICOM from sample: {out_path}")
            return out_path
        except Exception as e:
            logger.warning(f"Failed to re-tag sample CXR ({e}), creating synthetic")

    # Synthetic CXR fallback
    sop_uid = generate_uid()
    file_meta = pydicom.dataset.FileMetaDataset()
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.1.1"  # Digital X-Ray
    file_meta.MediaStorageSOPInstanceUID = sop_uid

    ds = pydicom.dataset.FileDataset(
        filename_or_obj="",
        dataset={},
        file_meta=file_meta,
        preamble=b"\x00" * 128,
    )

    ds.PatientName = dicom_patient_name
    ds.PatientID = patient_id
    ds.PatientBirthDate = MARIA_BIRTH_DATE
    ds.PatientSex = MARIA_SEX
    ds.PatientWeight = MARIA_WEIGHT
    ds.PatientAge = f"{MARIA_AGE:03d}Y"
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = generate_uid()
    ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.1.1"
    ds.SOPInstanceUID = sop_uid
    ds.StudyDate = STUDY_DATE
    ds.StudyTime = STUDY_TIME
    ds.StudyDescription = "CHEST PA AND LATERAL"
    ds.SeriesDescription = "PA UPRIGHT"
    ds.Modality = "CR"
    ds.BodyPartExamined = "CHEST"
    ds.InstitutionName = "Community Hospital - Rural NM"
    ds.Manufacturer = "HCLS AI Factory Demo"
    ds.InstanceNumber = 1
    ds.SeriesNumber = 1

    # Synthetic 512x512 CXR-like image
    rng = np.random.RandomState(99)
    cxr_img = rng.normal(2048, 200, (512, 512)).astype(np.uint16)
    cxr_img = np.clip(cxr_img, 0, 4095)

    ds.Rows = 512
    ds.Columns = 512
    ds.BitsAllocated = 16
    ds.BitsStored = 12
    ds.HighBit = 11
    ds.PixelRepresentation = 0
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.PixelData = cxr_img.astype(np.uint16).tobytes()

    out_path = cxr_dir / "cxr_001.dcm"
    ds.save_as(str(out_path), write_like_original=False)
    logger.info(f"Generated synthetic CXR DICOM: {out_path}")
    return out_path


# ═══════════════════════════════════════════════════════════════════════
# 3. PRE-COMPUTED WORKFLOW RESULTS
# ═══════════════════════════════════════════════════════════════════════


def precompute_workflow_results(output_dir: Path) -> Path:
    """Run all 9 workflows in mock mode, save results as JSON.

    Files:
        data/demo/results/ct_head_hemorrhage_result.json
        data/demo/results/ct_chest_lung_nodule_result.json
        data/demo/results/cxr_rapid_findings_result.json
        data/demo/results/ct_coronary_angiography_result.json
        data/demo/results/mri_brain_ms_lesion_result.json
        data/demo/results/mri_prostate_pirads_result.json
        data/demo/results/breast_birads_result.json
        data/demo/results/thyroid_tirads_result.json
        data/demo/results/liver_lirads_result.json

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the results directory.
    """
    from src.workflows import WORKFLOW_REGISTRY

    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    for wf_name, wf_class in WORKFLOW_REGISTRY.items():
        logger.info(f"Running workflow: {wf_name} (mock mode)")
        try:
            wf = wf_class(mock_mode=True)
            result = wf.run(input_path="")
            result_dict = result.model_dump(mode="json")
            result_dict["generated_at"] = datetime.now().isoformat()
            result_dict["demo_patient_id"] = DEFAULT_PATIENT_ID
            result_dict["demo_patient_name"] = DEFAULT_PATIENT_NAME

            out_path = results_dir / f"{wf_name}_result.json"
            out_path.write_text(json.dumps(result_dict, indent=2, default=str))
            logger.info(f"  Saved: {out_path.name} ({result.severity.value})")
        except Exception as e:
            logger.error(f"  Failed to run {wf_name}: {e}")

    return results_dir


# ═══════════════════════════════════════════════════════════════════════
# 4. PRE-COMPUTED RADIOMICS
# ═══════════════════════════════════════════════════════════════════════


def precompute_radiomics(output_dir: Path) -> Path:
    """Generate radiomics features in mock mode.

    Files:
        data/demo/radiomics/lung_nodule_features.json -- ~93 features for the nodule
        data/demo/radiomics/longitudinal_comparison.json -- comparison with prior
        data/demo/radiomics/feature_summary.txt -- embedding-ready text summary

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the radiomics directory.
    """
    radiomics_dir = output_dir / "radiomics"
    radiomics_dir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(42)

    # -- Lung nodule features (all 7 feature classes) --
    features: Dict[str, Any] = {
        "extraction_info": {
            "patient_id": DEFAULT_PATIENT_ID,
            "region": "right_upper_lobe_nodule",
            "modality": "CT",
            "voxel_spacing_mm": [0.703, 0.703, 1.25],
            "roi_shape_voxels": [12, 10, 8],
            "feature_classes": ["firstorder", "shape", "glcm", "glrlm", "glszm", "ngtdm", "gldm"],
            "extraction_mode": "mock",
            "generated_at": datetime.now().isoformat(),
        },
        "firstorder": {
            "Mean": round(rng.uniform(40.0, 80.0), 2),
            "Median": round(rng.uniform(38.0, 75.0), 2),
            "StandardDeviation": round(rng.uniform(12.0, 30.0), 2),
            "Variance": round(rng.uniform(144.0, 900.0), 2),
            "Skewness": round(rng.uniform(-0.5, 1.2), 4),
            "Kurtosis": round(rng.uniform(2.5, 5.0), 4),
            "Entropy": round(rng.uniform(4.0, 5.5), 4),
            "Energy": round(rng.uniform(1e7, 5e8), 0),
            "Minimum": round(rng.uniform(-100.0, 10.0), 1),
            "Maximum": round(rng.uniform(120.0, 250.0), 1),
            "Range": round(rng.uniform(200.0, 350.0), 1),
            "10Percentile": round(rng.uniform(15.0, 40.0), 1),
            "90Percentile": round(rng.uniform(90.0, 150.0), 1),
            "MeanAbsoluteDeviation": round(rng.uniform(8.0, 22.0), 2),
            "RobustMeanAbsoluteDeviation": round(rng.uniform(5.0, 18.0), 2),
            "RootMeanSquared": round(rng.uniform(45.0, 95.0), 2),
            "Uniformity": round(rng.uniform(0.03, 0.10), 4),
            "TotalEnergy": round(rng.uniform(5e7, 5e9), 0),
            "InterquartileRange": round(rng.uniform(20.0, 55.0), 1),
        },
        "shape": {
            "VoxelVolume": round(rng.uniform(800.0, 2500.0), 1),
            "MeshVolume": round(rng.uniform(780.0, 2450.0), 1),
            "SurfaceArea": round(rng.uniform(400.0, 1200.0), 1),
            "SurfaceVolumeRatio": round(rng.uniform(0.15, 0.45), 4),
            "Sphericity": round(rng.uniform(0.55, 0.85), 4),
            "Maximum3DDiameter": round(rng.uniform(8.0, 18.0), 1),
            "Maximum2DDiameterSlice": round(rng.uniform(6.0, 15.0), 1),
            "Maximum2DDiameterColumn": round(rng.uniform(5.0, 14.0), 1),
            "Maximum2DDiameterRow": round(rng.uniform(5.0, 13.0), 1),
            "MajorAxisLength": round(rng.uniform(8.0, 16.0), 2),
            "MinorAxisLength": round(rng.uniform(5.0, 12.0), 2),
            "LeastAxisLength": round(rng.uniform(4.0, 10.0), 2),
            "Elongation": round(rng.uniform(0.55, 0.90), 4),
            "Flatness": round(rng.uniform(0.40, 0.80), 4),
        },
        "glcm": {
            "Autocorrelation": round(rng.uniform(40.0, 200.0), 2),
            "ClusterProminence": round(rng.uniform(5000.0, 500000.0), 0),
            "ClusterShade": round(rng.uniform(-10000.0, 10000.0), 1),
            "ClusterTendency": round(rng.uniform(20.0, 200.0), 2),
            "Contrast": round(rng.uniform(8.0, 80.0), 2),
            "Correlation": round(rng.uniform(0.50, 0.92), 4),
            "DifferenceAverage": round(rng.uniform(3.0, 10.0), 3),
            "DifferenceEntropy": round(rng.uniform(2.0, 3.8), 4),
            "DifferenceVariance": round(rng.uniform(5.0, 40.0), 2),
            "JointAverage": round(rng.uniform(8.0, 35.0), 2),
            "JointEnergy": round(rng.uniform(0.005, 0.05), 5),
            "JointEntropy": round(rng.uniform(5.0, 7.5), 4),
            "MaximumProbability": round(rng.uniform(0.02, 0.10), 4),
            "SumAverage": round(rng.uniform(15.0, 70.0), 2),
            "SumEntropy": round(rng.uniform(3.5, 5.5), 4),
            "SumSquares": round(rng.uniform(10.0, 80.0), 2),
        },
        "glrlm": {
            "GrayLevelNonUniformity": round(rng.uniform(100.0, 2000.0), 1),
            "GrayLevelNonUniformityNormalized": round(rng.uniform(0.03, 0.10), 4),
            "HighGrayLevelRunEmphasis": round(rng.uniform(80.0, 800.0), 1),
            "LongRunEmphasis": round(rng.uniform(1.8, 3.5), 3),
            "RunEntropy": round(rng.uniform(3.5, 5.5), 4),
            "RunLengthNonUniformity": round(rng.uniform(800.0, 20000.0), 0),
            "RunPercentage": round(rng.uniform(0.75, 0.94), 4),
            "ShortRunEmphasis": round(rng.uniform(0.78, 0.95), 4),
        },
        "glszm": {
            "GrayLevelNonUniformity": round(rng.uniform(30.0, 500.0), 1),
            "HighGrayLevelZoneEmphasis": round(rng.uniform(80.0, 800.0), 1),
            "LargeAreaEmphasis": round(rng.uniform(200.0, 20000.0), 0),
            "SmallAreaEmphasis": round(rng.uniform(0.40, 0.85), 4),
            "ZoneEntropy": round(rng.uniform(5.0, 7.5), 4),
            "ZonePercentage": round(rng.uniform(0.01, 0.20), 4),
        },
        "ngtdm": {
            "Busyness": round(rng.uniform(0.5, 15.0), 3),
            "Coarseness": round(rng.uniform(0.0001, 0.01), 6),
            "Complexity": round(rng.uniform(5000.0, 100000.0), 0),
            "Contrast": round(rng.uniform(0.03, 0.20), 4),
            "Strength": round(rng.uniform(1.0, 10.0), 3),
        },
        "gldm": {
            "DependenceEntropy": round(rng.uniform(4.5, 6.5), 4),
            "DependenceNonUniformity": round(rng.uniform(200.0, 5000.0), 0),
            "GrayLevelNonUniformity": round(rng.uniform(30.0, 500.0), 1),
            "HighGrayLevelEmphasis": round(rng.uniform(80.0, 800.0), 1),
            "SmallDependenceEmphasis": round(rng.uniform(0.40, 0.85), 4),
        },
    }

    total_features = sum(len(v) for k, v in features.items() if isinstance(v, dict) and k != "extraction_info")
    features["extraction_info"]["total_features"] = total_features

    (radiomics_dir / "lung_nodule_features.json").write_text(
        json.dumps(features, indent=2)
    )
    logger.info(f"  Generated {total_features} radiomic features")

    # -- Longitudinal comparison --
    longitudinal = {
        "patient_id": DEFAULT_PATIENT_ID,
        "region": "right_upper_lobe_nodule",
        "current_study_date": TODAY.strftime("%Y-%m-%d"),
        "prior_study_date": (TODAY - timedelta(days=365)).strftime("%Y-%m-%d"),
        "interval_days": 365,
        "current": {
            "diameter_mm": 8.0,
            "volume_mm3": 268.1,
            "solid_component_mm": 4.0,
            "mean_hu": 42.5,
            "type": "part-solid",
        },
        "prior": {
            "diameter_mm": 5.0,
            "volume_mm3": 65.4,
            "solid_component_mm": 0.0,
            "mean_hu": -380.0,
            "type": "ground-glass",
        },
        "changes": {
            "diameter_change_mm": 3.0,
            "diameter_change_pct": 60.0,
            "volume_change_mm3": 202.7,
            "volume_change_pct": 309.9,
            "volume_doubling_time_days": 280,
            "morphology_change": "ground-glass to part-solid",
            "solid_component_emerged": True,
            "density_change_hu": 422.5,
        },
        "assessment": {
            "growth_rate": "concerning",
            "classification_change": "Lung-RADS 3 -> Lung-RADS 4B",
            "recommendation": "Tissue sampling recommended. Consider PET/CT.",
        },
    }
    (radiomics_dir / "longitudinal_comparison.json").write_text(
        json.dumps(longitudinal, indent=2)
    )

    # -- Embedding-ready text summary --
    summary_text = (
        f"Radiomic feature summary for {DEFAULT_PATIENT_ID} "
        f"right upper lobe nodule. "
        f"8mm part-solid nodule with 4mm solid component. "
        f"Volume 268.1 mm3, prior volume 65.4 mm3 (365 days ago). "
        f"Volume doubling time 280 days. "
        f"Morphology change: ground-glass to part-solid. "
        f"First-order mean HU 42.5, entropy 4.8. "
        f"Shape sphericity 0.72, elongation 0.68. "
        f"GLCM correlation 0.78, contrast 25.4. "
        f"Growth pattern and morphological evolution consistent with "
        f"adenocarcinoma spectrum. Lung-RADS 4B. "
        f"Tissue sampling recommended."
    )
    (radiomics_dir / "feature_summary.txt").write_text(summary_text)

    logger.info(f"  Generated radiomics data in {radiomics_dir}")
    return radiomics_dir


# ═══════════════════════════════════════════════════════════════════════
# 5. SAMPLE RADIOLOGY REPORTS
# ═══════════════════════════════════════════════════════════════════════


def generate_sample_reports(output_dir: Path) -> Path:
    """Generate sample radiology reports for demo.

    Files:
        data/demo/reports/ct_chest_maria_santos.txt -- the main demo report
        data/demo/reports/ct_chest_prior.txt -- prior study report for comparison
        data/demo/reports/cxr_pneumonia.txt -- sample CXR report

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the reports directory.
    """
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Main CT chest report (Maria Santos)
    prior_date = (TODAY - timedelta(days=365)).strftime("%Y-%m-%d")
    ct_report = f"""\
EXAMINATION: CT Chest with Contrast

CLINICAL HISTORY: 58-year-old female, 20 pack-year smoking history, annual lung cancer screening. Prior CT Chest {prior_date} showed 5mm ground-glass nodule right upper lobe.

TECHNIQUE: Helical CT, 1.25mm slice thickness, 80mL Omnipaque 350 IV contrast, arterial phase. kVp 120, effective mAs 200.

COMPARISON: CT Chest dated {prior_date}.

FINDINGS:
Lungs: An 8mm part-solid nodule is present in the right upper lobe (series 3, image 127), previously measuring 5mm on the prior study dated {prior_date}. The solid component now measures 4mm. Volume doubling time is approximately 280 days. No new pulmonary nodules are identified. No consolidation, ground-glass opacity, or air trapping.
Airways: The trachea and main bronchi are patent without evidence of endobronchial lesion.
Mediastinum: No pathologically enlarged mediastinal or hilar lymph nodes. The largest prevascular lymph node measures 8mm in short axis, stable. Heart size is normal. No pericardial effusion.
Pleura: No pleural effusion. No pneumothorax.
Upper Abdomen: Limited evaluation. The visualized liver, spleen, and adrenal glands appear unremarkable.
Bones: No suspicious osseous lesions. Mild degenerative changes of the thoracic spine.

IMPRESSION:
1. Growing 8mm part-solid right upper lobe nodule with 4mm solid component, previously 5mm. Volume doubling time 280 days. Lung-RADS 4B. Recommend tissue sampling. Consider PET/CT for further characterization.
2. Stable 8mm prevascular lymph node.
3. Otherwise unremarkable CT chest.

Electronically signed by:
Elena Rodriguez, MD
Board Certified Radiologist
Community Hospital -- Rural NM
"""

    (reports_dir / "ct_chest_maria_santos.txt").write_text(ct_report)

    # Prior study report
    prior_report_date = (TODAY - timedelta(days=365)).strftime("%Y-%m-%d")
    two_years_ago = (TODAY - timedelta(days=730)).strftime("%Y-%m-%d")
    prior_report = f"""\
EXAMINATION: CT Chest Low-Dose Lung Cancer Screening

CLINICAL HISTORY: 57-year-old female, 20 pack-year smoking history. Annual lung cancer screening per USPSTF guidelines. No prior CT chest available for comparison.

TECHNIQUE: Low-dose helical CT, 1.25mm slice thickness, no IV contrast. kVp 120, effective mAs 30.

COMPARISON: None.

FINDINGS:
Lungs: A 5mm ground-glass nodule is present in the right upper lobe, posterior segment (series 2, image 98). No solid component identified. No other pulmonary nodules. No consolidation or air trapping. Mild centrilobular emphysema in bilateral upper lobes.
Airways: Patent trachea and main bronchi. No endobronchial lesion.
Mediastinum: No enlarged mediastinal or hilar lymph nodes. An 8mm prevascular lymph node is noted, likely reactive. Normal cardiac silhouette.
Pleura: No effusion. No pneumothorax.
Upper Abdomen: Limited evaluation. Unremarkable.
Bones: Mild thoracic spondylosis. No suspicious osseous lesions.

IMPRESSION:
1. 5mm ground-glass nodule, right upper lobe. Lung-RADS 3 -- probably benign. Recommend 6-month follow-up low-dose CT.
2. Mild centrilobular emphysema, bilateral upper lobes.
3. 8mm prevascular lymph node, likely reactive.

Electronically signed by:
James Chen, MD
Board Certified Radiologist
Community Hospital -- Rural NM
"""

    (reports_dir / "ct_chest_prior.txt").write_text(prior_report)

    # CXR pneumonia report
    cxr_report = """\
EXAMINATION: Portable AP Chest X-Ray

CLINICAL HISTORY: 45-year-old female with 5 days of productive cough, fever 39.2C, dyspnea, pleuritic chest pain. SpO2 91% on room air. WBC 18,200. Procalcitonin 2.4 ng/mL. Rule out pneumonia.

TECHNIQUE: Single AP portable radiograph.

COMPARISON: None available.

FINDINGS:
Lungs: Dense lobar consolidation with air bronchograms in the right lower lobe. Patchy consolidation with ground-glass opacity in the left lower lobe. No pneumothorax. No upper lobe infiltrates.
Pleura: Small right-sided pleural effusion with meniscus sign at the right costophrenic angle. Left costophrenic angle is clear.
Mediastinum: Normal cardiac silhouette, cardiothoracic ratio 0.47. Mediastinal contours are within normal limits.
Bones: No acute osseous abnormality.

IMPRESSION:
1. Bilateral lower lobe pneumonia, more severe on the right with dense lobar consolidation and air bronchograms. Small right pleural effusion.
2. Clinical correlation with elevated WBC and procalcitonin suggests bacterial etiology. Consider blood cultures and empiric antibiotics per institutional protocol.
3. Recommend follow-up imaging in 4-6 weeks to document resolution.

Electronically signed by:
Sarah Kim, MD
Board Certified Radiologist
Community Hospital -- Rural NM
"""

    (reports_dir / "cxr_pneumonia.txt").write_text(cxr_report)

    logger.info(f"  Generated 3 radiology reports in {reports_dir}")
    return reports_dir


# ═══════════════════════════════════════════════════════════════════════
# 6. PARSED REPORT NLP RESULTS
# ═══════════════════════════════════════════════════════════════════════


def precompute_report_nlp(output_dir: Path) -> Path:
    """Parse the sample reports and save results.

    Files:
        data/demo/parsed_reports/ct_chest_parsed.json -- sections, entities, measurements

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the parsed reports directory.
    """
    parsed_dir = output_dir / "parsed_reports"
    parsed_dir.mkdir(parents=True, exist_ok=True)

    prior_date = (TODAY - timedelta(days=365)).strftime("%Y-%m-%d")

    parsed = {
        "patient_id": DEFAULT_PATIENT_ID,
        "patient_name": DEFAULT_PATIENT_NAME,
        "report_type": "CT Chest with Contrast",
        "parsed_at": datetime.now().isoformat(),
        "sections": {
            "examination": "CT Chest with Contrast",
            "clinical_history": (
                f"58-year-old female, 20 pack-year smoking history, annual lung cancer "
                f"screening. Prior CT Chest {prior_date} showed 5mm ground-glass "
                f"nodule right upper lobe."
            ),
            "technique": (
                "Helical CT, 1.25mm slice thickness, 80mL Omnipaque 350 IV contrast, "
                "arterial phase. kVp 120, effective mAs 200."
            ),
            "comparison": f"CT Chest dated {prior_date}.",
            "findings": (
                "Lungs: An 8mm part-solid nodule is present in the right upper lobe "
                "(series 3, image 127), previously measuring 5mm on the prior study. "
                "The solid component now measures 4mm. Volume doubling time is "
                "approximately 280 days. No new pulmonary nodules. "
                "Airways: Patent. Mediastinum: No enlarged lymph nodes. 8mm prevascular "
                "lymph node, stable. Normal heart size. "
                "Pleura: No effusion. No pneumothorax. "
                "Upper Abdomen: Unremarkable. "
                "Bones: Mild degenerative changes thoracic spine."
            ),
            "impression": (
                "1. Growing 8mm part-solid RUL nodule with 4mm solid component, "
                "previously 5mm. VDT 280 days. Lung-RADS 4B. Recommend tissue "
                "sampling. Consider PET/CT. "
                "2. Stable 8mm prevascular lymph node. "
                "3. Otherwise unremarkable CT chest."
            ),
        },
        "entities": [
            {"entity_type": "finding", "value": "part-solid nodule", "context": "8mm part-solid nodule in right upper lobe"},
            {"entity_type": "anatomy", "value": "right upper lobe", "context": "nodule in the right upper lobe"},
            {"entity_type": "finding", "value": "solid component", "context": "solid component now measures 4mm"},
            {"entity_type": "finding", "value": "prevascular lymph node", "context": "8mm prevascular lymph node, stable"},
            {"entity_type": "finding", "value": "degenerative changes", "context": "mild degenerative changes thoracic spine"},
            {"entity_type": "anatomy", "value": "thoracic spine", "context": "degenerative changes of the thoracic spine"},
            {"entity_type": "laterality", "value": "right", "context": "right upper lobe"},
            {"entity_type": "recommendation", "value": "tissue sampling", "context": "Recommend tissue sampling"},
            {"entity_type": "recommendation", "value": "PET/CT", "context": "Consider PET/CT for further characterization"},
            {"entity_type": "classification", "value": "Lung-RADS 4B", "context": "Lung-RADS 4B"},
        ],
        "measurements": [
            {"value": 8.0, "unit": "mm", "finding": "nodule", "qualifier": "part-solid", "location": "right upper lobe"},
            {"value": 5.0, "unit": "mm", "finding": "nodule", "qualifier": "ground-glass", "location": "right upper lobe (prior)"},
            {"value": 4.0, "unit": "mm", "finding": "solid component", "qualifier": None, "location": "right upper lobe"},
            {"value": 280.0, "unit": "days", "finding": "volume doubling time", "qualifier": None, "location": None},
            {"value": 8.0, "unit": "mm", "finding": "lymph node", "qualifier": "prevascular", "location": "mediastinum"},
            {"value": 1.25, "unit": "mm", "finding": "slice thickness", "qualifier": None, "location": None},
            {"value": 120.0, "unit": "kVp", "finding": "tube voltage", "qualifier": None, "location": None},
        ],
        "critical_finding": True,
        "modality": "CT",
        "body_region": "chest",
        "comparison_date": prior_date,
    }

    (parsed_dir / "ct_chest_parsed.json").write_text(json.dumps(parsed, indent=2))
    logger.info(f"  Generated parsed report NLP results in {parsed_dir}")
    return parsed_dir


# ═══════════════════════════════════════════════════════════════════════
# 7. CROSS-MODAL TRIGGER RESULTS
# ═══════════════════════════════════════════════════════════════════════


def precompute_cross_modal(output_dir: Path) -> Path:
    """Generate cross-modal trigger results for the lung nodule case.

    Files:
        data/demo/cross_modal/lung_rads_4b_trigger.json -- trigger activation
        data/demo/cross_modal/egfr_evidence.json -- simulated genomic evidence

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the cross-modal directory.
    """
    cross_modal_dir = output_dir / "cross_modal"
    cross_modal_dir.mkdir(parents=True, exist_ok=True)

    # Trigger activation
    trigger = {
        "trigger_id": "TRIG-001",
        "trigger_type": "lung_rads_threshold",
        "source_workflow": "ct_chest_lung_nodule",
        "source_classification": "Lung-RADS 4B",
        "source_severity": "critical",
        "activation_rule": "Lung-RADS >= 4A triggers lung cancer genomic query",
        "activated_at": datetime.now().isoformat(),
        "patient_id": DEFAULT_PATIENT_ID,
        "finding_summary": "8mm part-solid nodule, right upper lobe, VDT 280 days",
        "genomic_queries": [
            "lung cancer driver mutations EGFR ALK ROS1 KRAS",
            "non-small cell lung cancer NSCLC targeted therapy genomics",
            "lung adenocarcinoma molecular subtypes precision medicine",
        ],
        "target_genes": ["EGFR", "ALK", "ROS1", "KRAS", "BRAF", "MET"],
        "enrichment_summary": (
            "Cross-modal genomic enrichment identified actionable driver "
            "mutations relevant to lung adenocarcinoma. EGFR L858R (exon 21 "
            "missense) found with ClinVar pathogenic classification and high "
            "AlphaMissense score (0.94). No actionable ALK, ROS1, or KRAS "
            "variants detected. Recommend molecular-guided therapy planning."
        ),
        "genomic_hit_count": 12,
        "query_count": 3,
    }
    (cross_modal_dir / "lung_rads_4b_trigger.json").write_text(
        json.dumps(trigger, indent=2)
    )

    # EGFR evidence
    egfr_evidence = {
        "patient_id": DEFAULT_PATIENT_ID,
        "trigger_id": "TRIG-001",
        "queried_at": datetime.now().isoformat(),
        "genes_queried": ["EGFR", "ALK", "ROS1", "KRAS"],
        "results": [
            {
                "gene": "EGFR",
                "variant": "L858R",
                "hgvs_p": "p.Leu858Arg",
                "hgvs_c": "c.2573T>G",
                "exon": 21,
                "variant_type": "missense",
                "clinvar_significance": "Pathogenic",
                "clinvar_review_status": "reviewed by expert panel",
                "clinvar_id": "VCV000016609",
                "alphamissense_score": 0.94,
                "alphamissense_class": "likely_pathogenic",
                "cosmic_id": "COSM6224",
                "cosmic_frequency": 0.29,
                "population_frequency_gnomad": 0.0,
                "functional_impact": "Constitutive activation of EGFR tyrosine kinase domain",
                "therapeutic_relevance": [
                    "Sensitive to EGFR TKIs: erlotinib, gefitinib, afatinib",
                    "First-line osimertinib per NCCN NSCLC guidelines",
                    "Median PFS 18.9 months with osimertinib (FLAURA trial)",
                ],
                "clinical_trials": [
                    "NCT02296125 (FLAURA)",
                    "NCT04487080 (LAURA)",
                    "NCT03944772 (ADAURA)",
                ],
                "actionable": True,
                "evidence_level": "1A",
            },
            {
                "gene": "ALK",
                "variant": None,
                "actionable": False,
                "summary": "No actionable ALK rearrangements or mutations found",
            },
            {
                "gene": "ROS1",
                "variant": None,
                "actionable": False,
                "summary": "No actionable ROS1 rearrangements or mutations found",
            },
            {
                "gene": "KRAS",
                "variant": None,
                "actionable": False,
                "summary": "No actionable KRAS mutations found (G12C negative)",
            },
        ],
        "summary": (
            "Genomic analysis reveals EGFR L858R activating mutation (exon 21 "
            "missense, ClinVar pathogenic, AlphaMissense 0.94). This is the "
            "most common EGFR sensitizing mutation in NSCLC, conferring "
            "sensitivity to EGFR tyrosine kinase inhibitors. First-line "
            "osimertinib recommended per NCCN guidelines. No ALK, ROS1, or "
            "KRAS actionable variants detected."
        ),
    }
    (cross_modal_dir / "egfr_evidence.json").write_text(
        json.dumps(egfr_evidence, indent=2)
    )

    logger.info(f"  Generated cross-modal trigger data in {cross_modal_dir}")
    return cross_modal_dir


# ═══════════════════════════════════════════════════════════════════════
# 8. ANALYTICS DEMO DATA
# ═══════════════════════════════════════════════════════════════════════


def precompute_analytics(output_dir: Path) -> Path:
    """Generate analytics demo data.

    Files:
        data/demo/analytics/population_500_studies.json -- 500 synthetic studies
        data/demo/analytics/population_summary.json -- pre-computed summary

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the analytics directory.
    """
    analytics_dir = output_dir / "analytics"
    analytics_dir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(12345)

    modalities = ["CT", "MRI", "CXR", "US", "PET-CT", "Mammography"]
    modality_weights = [0.35, 0.20, 0.25, 0.10, 0.05, 0.05]

    body_regions = ["chest", "head", "abdomen", "brain", "breast", "cardiac", "pelvis", "spine"]
    region_weights = [0.25, 0.15, 0.15, 0.10, 0.08, 0.10, 0.10, 0.07]

    severities = ["normal", "routine", "significant", "urgent", "critical"]
    severity_weights = [0.30, 0.25, 0.20, 0.15, 0.10]

    finding_types = [
        "nodule", "consolidation", "effusion", "hemorrhage", "fracture",
        "mass", "calcification", "stenosis", "lesion", "normal",
        "atelectasis", "edema", "pneumothorax", "lymphadenopathy",
    ]

    workflows = [
        "ct_head_hemorrhage", "ct_chest_lung_nodule", "cxr_rapid_findings",
        "ct_coronary_angiography", "mri_brain_ms_lesion", "mri_prostate_pirads",
        "breast_birads", "thyroid_tirads", "liver_lirads",
    ]

    studies = []
    start_date = TODAY - timedelta(days=90)

    for i in range(500):
        study_dt = start_date + timedelta(
            days=rng.randint(0, 90),
            hours=rng.randint(6, 22),
            minutes=rng.randint(0, 59),
        )
        modality = rng.choices(modalities, weights=modality_weights, k=1)[0]
        region = rng.choices(body_regions, weights=region_weights, k=1)[0]
        severity = rng.choices(severities, weights=severity_weights, k=1)[0]
        proc_time = round(rng.gauss(1800, 600), 0)
        proc_time = max(200, min(8000, proc_time))

        age = rng.randint(25, 85)
        sex = rng.choice(["M", "F"])
        finding = rng.choice(finding_types)

        studies.append({
            "study_id": f"STU-{i + 1:04d}",
            "patient_id": f"PAT-{rng.randint(1000, 9999)}",
            "study_date": study_dt.strftime("%Y-%m-%d"),
            "study_time": study_dt.strftime("%H:%M:%S"),
            "modality": modality,
            "body_region": region,
            "workflow": rng.choice(workflows),
            "severity": severity,
            "primary_finding": finding,
            "processing_time_ms": proc_time,
            "ai_confidence": round(rng.uniform(0.65, 0.99), 3),
            "patient_age": age,
            "patient_sex": sex,
            "critical_finding_alert": severity in ("critical", "urgent"),
        })

    (analytics_dir / "population_500_studies.json").write_text(
        json.dumps(studies, indent=2)
    )

    # Pre-computed summary
    modality_dist = {}
    region_dist = {}
    severity_dist = {}
    finding_counts: Dict[str, int] = {}
    total_proc_time = 0.0
    critical_count = 0

    for s in studies:
        modality_dist[s["modality"]] = modality_dist.get(s["modality"], 0) + 1
        region_dist[s["body_region"]] = region_dist.get(s["body_region"], 0) + 1
        severity_dist[s["severity"]] = severity_dist.get(s["severity"], 0) + 1
        finding_counts[s["primary_finding"]] = finding_counts.get(s["primary_finding"], 0) + 1
        total_proc_time += s["processing_time_ms"]
        if s["severity"] in ("critical", "urgent"):
            critical_count += 1

    finding_prevalence = {k: round(v / 500 * 100, 1) for k, v in finding_counts.items()}

    summary = {
        "total_studies": 500,
        "modality_distribution": modality_dist,
        "body_region_distribution": region_dist,
        "severity_distribution": severity_dist,
        "finding_prevalence": finding_prevalence,
        "mean_processing_time_ms": round(total_proc_time / 500, 1),
        "date_range": {
            "start": start_date.strftime("%Y-%m-%d"),
            "end": TODAY.strftime("%Y-%m-%d"),
        },
        "studies_with_critical_findings": critical_count,
        "critical_finding_rate": round(critical_count / 500 * 100, 1),
        "generated_at": datetime.now().isoformat(),
    }

    (analytics_dir / "population_summary.json").write_text(
        json.dumps(summary, indent=2)
    )

    logger.info(f"  Generated 500 synthetic studies + summary in {analytics_dir}")
    return analytics_dir


# ═══════════════════════════════════════════════════════════════════════
# 9. DICOM SR
# ═══════════════════════════════════════════════════════════════════════


def precompute_dicom_sr(output_dir: Path) -> Path:
    """Generate DICOM SR from the lung nodule workflow result.

    Attempts to use the DICOMSRExporter if highdicom is available.
    Falls back to a metadata-only JSON if not.

    Files:
        data/demo/dicom_sr/ct_chest_lung_nodule_sr.dcm -- TID 1500 Measurement Report
        (or .json fallback if highdicom is unavailable)

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the DICOM SR directory.
    """
    sr_dir = output_dir / "dicom_sr"
    sr_dir.mkdir(parents=True, exist_ok=True)

    # Load the pre-computed workflow result
    result_path = output_dir / "results" / "ct_chest_lung_nodule_result.json"
    if not result_path.exists():
        logger.warning("Lung nodule result not found; generating workflow first")
        precompute_workflow_results(output_dir)

    from src.models import WorkflowResult

    result_data = json.loads(result_path.read_text())
    wf_result = WorkflowResult(**{
        k: v for k, v in result_data.items()
        if k in WorkflowResult.model_fields
    })

    # Try highdicom-based export
    try:
        from src.export_dicom_sr import DICOMSRExporter, DICOMSRConfig

        config = DICOMSRConfig(
            institution_name="Community Hospital - Rural NM",
            device_name="Clinical Imaging Engine v2.0",
            manufacturer="HCLS AI Factory (Apache 2.0)",
        )
        exporter = DICOMSRExporter(config=config)

        if exporter.available:
            sr_path = str(sr_dir / "ct_chest_lung_nodule_sr.dcm")
            sr_result = exporter.export_workflow_result(
                workflow_result=wf_result,
                output_path=sr_path,
                patient_id=DEFAULT_PATIENT_ID,
                patient_name=f"SANTOS^MARIA",
                study_date=STUDY_DATE,
                accession_number=f"ACC{DEFAULT_PATIENT_ID}",
            )
            if sr_result.success:
                logger.info(f"  Generated DICOM SR: {sr_path}")
                return sr_dir
            else:
                logger.warning(f"  SR export failed: {sr_result.error}")
        else:
            logger.info("  highdicom not available, generating JSON fallback")
    except Exception as e:
        logger.warning(f"  DICOM SR export failed ({e}), generating JSON fallback")

    # JSON fallback
    sr_json = {
        "type": "DICOM SR TID 1500 Measurement Report (JSON fallback)",
        "note": "Install highdicom for real DICOM SR generation",
        "patient_id": DEFAULT_PATIENT_ID,
        "patient_name": "SANTOS^MARIA",
        "study_date": STUDY_DATE,
        "workflow": "ct_chest_lung_nodule",
        "classification": wf_result.classification,
        "severity": wf_result.severity.value,
        "findings": wf_result.findings,
        "measurements": wf_result.measurements,
        "generated_at": datetime.now().isoformat(),
    }
    fallback_path = sr_dir / "ct_chest_lung_nodule_sr.json"
    fallback_path.write_text(json.dumps(sr_json, indent=2, default=str))
    logger.info(f"  Generated DICOM SR JSON fallback: {fallback_path}")
    return sr_dir


# ═══════════════════════════════════════════════════════════════════════
# 10. ENGINE 3 HANDOFF (TARGET HYPOTHESIS)
# ═══════════════════════════════════════════════════════════════════════


def generate_engine3_target(output_dir: Path) -> Path:
    """Generate target hypothesis for Engine 3 drug discovery handoff.

    Files:
        data/demo/engine3_handoff/egfr_target_hypothesis.json

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the engine3 handoff directory.
    """
    handoff_dir = output_dir / "engine3_handoff"
    handoff_dir.mkdir(parents=True, exist_ok=True)

    target = {
        "gene": "EGFR",
        "protein": "Epidermal Growth Factor Receptor",
        "uniprot_id": "P00533",
        "variant": "L858R",
        "variant_type": "missense",
        "exon": 21,
        "therapeutic_area": "Non-small cell lung cancer",
        "mechanism": "Kinase domain activating mutation",
        "confidence": "HIGH",
        "priority": 5,
        "pdb_ids": ["1M17", "4ZAU", "5CAL"],
        "reference_drug": "Erlotinib",
        "reference_smiles": "C=Cc1cccc(NC2=NC=NC3=CC(OCCOC)=C(OCCOC)C=C23)c1",
        "druggability": "high",
        "source": "Clinical Imaging Engine \u2014 Cross-modal trigger (Lung-RADS 4B)",
        "imaging_trigger": {
            "workflow": "ct_chest_lung_nodule",
            "classification": "Lung-RADS 4B",
            "finding": "8mm part-solid nodule, RUL",
            "cross_modal_genes": ["EGFR", "ALK", "ROS1", "KRAS"],
        },
        "generated_at": datetime.now().isoformat(),
        "pipeline_version": "2.0.0",
    }

    (handoff_dir / "egfr_target_hypothesis.json").write_text(
        json.dumps(target, indent=2)
    )
    logger.info(f"  Generated Engine 3 target hypothesis in {handoff_dir}")
    return handoff_dir


# ═══════════════════════════════════════════════════════════════════════
# 11. DOSE RECORDS
# ═══════════════════════════════════════════════════════════════════════


def generate_dose_records(output_dir: Path) -> Path:
    """Generate dose tracking records for Maria Santos.

    Includes 4 studies showing cumulative dose tracking across prior and
    current imaging studies.

    Files:
        data/demo/dose/maria_santos_dose_history.json

    Args:
        output_dir: Root output directory.

    Returns:
        Path to the dose directory.
    """
    dose_dir = output_dir / "dose"
    dose_dir.mkdir(parents=True, exist_ok=True)

    dose_history = {
        "patient_id": DEFAULT_PATIENT_ID,
        "patient_name": DEFAULT_PATIENT_NAME,
        "generated_at": datetime.now().isoformat(),
        "studies": [
            {
                "study_date": (TODAY - timedelta(days=730)).strftime("%Y-%m-%d"),
                "modality": "CT",
                "protocol": "CT Chest Low-Dose Screening",
                "body_region": "chest",
                "effective_dose_msv": 1.5,
                "dlp_mgy_cm": 75.0,
                "ctdi_vol_mgy": 2.5,
                "num_exposures": 1,
                "scanner_model": "GE Revolution CT",
                "indication": "Baseline lung cancer screening",
            },
            {
                "study_date": (TODAY - timedelta(days=365)).strftime("%Y-%m-%d"),
                "modality": "CT",
                "protocol": "CT Chest Low-Dose Screening",
                "body_region": "chest",
                "effective_dose_msv": 1.6,
                "dlp_mgy_cm": 80.0,
                "ctdi_vol_mgy": 2.7,
                "num_exposures": 1,
                "scanner_model": "GE Revolution CT",
                "indication": "Annual lung cancer screening (5mm GGN follow-up)",
            },
            {
                "study_date": (TODAY - timedelta(days=14)).strftime("%Y-%m-%d"),
                "modality": "XR",
                "protocol": "Chest PA and Lateral",
                "body_region": "chest",
                "effective_dose_msv": 0.06,
                "dap_mgy_cm2": 12.0,
                "num_exposures": 2,
                "scanner_model": "Philips DigitalDiagnost",
                "indication": "Pre-CT chest X-ray, cough evaluation",
            },
            {
                "study_date": TODAY.strftime("%Y-%m-%d"),
                "modality": "CT",
                "protocol": "CT Chest with Contrast",
                "body_region": "chest",
                "effective_dose_msv": 7.0,
                "dlp_mgy_cm": 420.0,
                "ctdi_vol_mgy": 12.5,
                "num_exposures": 1,
                "scanner_model": "GE Revolution CT",
                "indication": "Growing lung nodule evaluation, contrast-enhanced",
            },
        ],
        "cumulative": {
            "total_effective_dose_msv": 10.16,
            "study_count": 4,
            "date_range": {
                "first": (TODAY - timedelta(days=730)).strftime("%Y-%m-%d"),
                "last": TODAY.strftime("%Y-%m-%d"),
            },
            "by_modality": {
                "CT": 10.1,
                "XR": 0.06,
            },
            "by_body_region": {
                "chest": 10.16,
            },
            "alert_level": "normal",
            "alert_message": (
                "Cumulative dose 10.16 mSv over 2 years is within normal range "
                "for diagnostic imaging. Current study (7.0 mSv) is within DRL "
                "for CT Chest with Contrast (national DRL: 10 mSv). "
                "No dose optimization alert required."
            ),
        },
        "drl_comparison": {
            "current_study": {
                "protocol": "CT Chest with Contrast",
                "effective_dose_msv": 7.0,
                "national_drl_msv": 10.0,
                "pct_of_drl": 70.0,
                "status": "within_drl",
                "optimization_note": (
                    "Current dose is 70% of national DRL. Protocol is "
                    "appropriately optimized for diagnostic quality."
                ),
            },
        },
    }

    (dose_dir / "maria_santos_dose_history.json").write_text(
        json.dumps(dose_history, indent=2)
    )
    logger.info(f"  Generated dose records in {dose_dir}")
    return dose_dir


# ═══════════════════════════════════════════════════════════════════════
# MAIN ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description="Prepare Clinical Imaging Engine demo data"
    )
    parser.add_argument(
        "--output-dir", default="data/demo",
        help="Output directory (default: data/demo)",
    )
    parser.add_argument(
        "--patient-name", default=DEFAULT_PATIENT_NAME,
        help="Demo patient name (default: Maria Santos)",
    )
    parser.add_argument(
        "--patient-id", default=DEFAULT_PATIENT_ID,
        help="Demo patient ID (default: MS2026)",
    )
    parser.add_argument(
        "--num-ct-slices", type=int, default=DEFAULT_NUM_CT_SLICES,
        help="Number of CT slices to generate (default: 64)",
    )
    args = parser.parse_args()

    # Resolve output directory relative to project root
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = PROJECT_ROOT / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Clinical Imaging Engine \u2014 Demo Data Preparation")
    print("=" * 60)
    print(f"  Patient:    {args.patient_name} ({args.patient_id})")
    print(f"  Output:     {output_dir}")
    print(f"  CT slices:  {args.num_ct_slices}")
    print(f"  Date:       {TODAY.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    steps = [
        ("1/11", "DICOM CT Chest Study", lambda: generate_ct_chest_dicom_study(
            output_dir, args.patient_name, args.patient_id, args.num_ct_slices,
        )),
        ("2/11", "CXR DICOM Study", lambda: generate_cxr_dicom(
            output_dir, args.patient_name, args.patient_id,
        )),
        ("3/11", "Workflow Results (9 workflows)", lambda: precompute_workflow_results(output_dir)),
        ("4/11", "Radiomics Features", lambda: precompute_radiomics(output_dir)),
        ("5/11", "Sample Radiology Reports", lambda: generate_sample_reports(output_dir)),
        ("6/11", "Report NLP Results", lambda: precompute_report_nlp(output_dir)),
        ("7/11", "Cross-Modal Trigger Results", lambda: precompute_cross_modal(output_dir)),
        ("8/11", "Analytics Demo Data (500 studies)", lambda: precompute_analytics(output_dir)),
        ("9/11", "DICOM SR", lambda: precompute_dicom_sr(output_dir)),
        ("10/11", "Engine 3 Target Hypothesis", lambda: generate_engine3_target(output_dir)),
        ("11/11", "Dose Records", lambda: generate_dose_records(output_dir)),
    ]

    results = []
    total_start = time.time()

    for step_num, step_name, step_fn in steps:
        print(f"\n[{step_num}] {step_name}...")
        step_start = time.time()
        try:
            result_path = step_fn()
            elapsed = time.time() - step_start
            results.append((step_name, "OK", elapsed, result_path))
            print(f"  Done ({elapsed:.1f}s)")
        except Exception as e:
            elapsed = time.time() - step_start
            results.append((step_name, f"FAILED: {e}", elapsed, None))
            logger.error(f"  Failed: {e}")
            print(f"  FAILED ({elapsed:.1f}s): {e}")

    total_elapsed = time.time() - total_start

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    ok_count = sum(1 for _, status, _, _ in results if status == "OK")
    fail_count = len(results) - ok_count

    for step_name, status, elapsed, result_path in results:
        icon = "OK" if status == "OK" else "FAIL"
        print(f"  [{icon}] {step_name} ({elapsed:.1f}s)")

    # Count generated files
    total_files = 0
    total_size = 0
    if output_dir.exists():
        for f in output_dir.rglob("*"):
            if f.is_file():
                total_files += 1
                total_size += f.stat().st_size

    print(f"\n  Steps:      {ok_count} OK, {fail_count} failed")
    print(f"  Files:      {total_files}")
    print(f"  Total size: {total_size / (1024 * 1024):.1f} MB")
    print(f"  Time:       {total_elapsed:.1f}s")
    print(f"  Output:     {output_dir}")
    print("=" * 60)

    if fail_count > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
