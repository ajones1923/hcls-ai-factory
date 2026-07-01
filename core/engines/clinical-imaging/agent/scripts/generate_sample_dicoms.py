#!/usr/bin/env python3
"""Generate synthetic DICOM series from existing NIfTI sample volumes.

Creates properly formatted DICOM series with valid UIDs and metadata,
suitable for upload to Orthanc for OHIF viewer demos.

Usage:
    python scripts/generate_sample_dicoms.py                    # Generate only
    python scripts/generate_sample_dicoms.py --upload            # Generate + upload to Orthanc
    python scripts/generate_sample_dicoms.py --orthanc-url http://localhost:8042

Author: Adam Jones
Date: March 2026
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def generate_dicom_series(
    nifti_path: Path,
    output_dir: Path,
    patient_name: str = "DEMO^Patient",
    patient_id: str = "DEMO-001",
    study_description: str = "Demo Study",
    series_description: str = "Demo Series",
    modality: str = "CT",
    max_slices: int = 50,
) -> Path:
    """Convert a NIfTI volume to a DICOM series.

    Args:
        nifti_path: Path to input NIfTI file (.nii or .nii.gz).
        output_dir: Directory to write DICOM files.
        patient_name: DICOM patient name tag.
        patient_id: DICOM patient ID tag.
        study_description: DICOM study description.
        series_description: DICOM series description.
        modality: Imaging modality (CT, MR, etc.).
        max_slices: Maximum number of slices to generate.

    Returns:
        Path to the output directory containing DICOM files.
    """
    import nibabel as nib
    import pydicom
    from pydicom.uid import generate_uid, ExplicitVRLittleEndian

    # Load NIfTI volume
    print(f"Loading {nifti_path}...")
    img = nib.load(str(nifti_path))
    data = img.get_fdata()

    # Take axial slices (along z-axis)
    n_slices = min(data.shape[2], max_slices)
    step = max(1, data.shape[2] // n_slices)
    slice_indices = list(range(0, data.shape[2], step))[:n_slices]

    # Normalize to 12-bit range (0-4095)
    data_min, data_max = float(np.nanmin(data)), float(np.nanmax(data))
    if data_max > data_min:
        normalized = ((data - data_min) / (data_max - data_min) * 4095).astype(np.uint16)
    else:
        normalized = np.zeros_like(data, dtype=np.uint16)

    # Generate UIDs
    study_uid = generate_uid()
    series_uid = generate_uid()
    frame_of_ref_uid = generate_uid()

    series_dir = output_dir / f"{patient_id}_{series_description.replace(' ', '_')}"
    series_dir.mkdir(parents=True, exist_ok=True)

    now = datetime.now()
    study_date = now.strftime("%Y%m%d")
    study_time = now.strftime("%H%M%S")

    print(f"Generating {len(slice_indices)} DICOM slices...")

    for i, z_idx in enumerate(slice_indices):
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

        # Patient
        ds.PatientName = patient_name
        ds.PatientID = patient_id
        ds.PatientBirthDate = "19700101"
        ds.PatientSex = "O"

        # Study
        ds.StudyInstanceUID = study_uid
        ds.StudyDate = study_date
        ds.StudyTime = study_time
        ds.StudyDescription = study_description
        ds.StudyID = "1"
        ds.AccessionNumber = ""

        # Series
        ds.SeriesInstanceUID = series_uid
        ds.SeriesNumber = 1
        ds.SeriesDescription = series_description
        ds.Modality = modality
        ds.FrameOfReferenceUID = frame_of_ref_uid

        # Instance
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        ds.SOPInstanceUID = sop_uid
        ds.InstanceNumber = i + 1
        ds.ImagePositionPatient = [0.0, 0.0, float(z_idx)]
        ds.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        ds.SliceLocation = float(z_idx)
        ds.SliceThickness = float(step)
        ds.PixelSpacing = [1.0, 1.0]

        # Pixel data
        slice_data = normalized[:, :, z_idx]
        ds.Rows = slice_data.shape[0]
        ds.Columns = slice_data.shape[1]
        ds.BitsAllocated = 16
        ds.BitsStored = 12
        ds.HighBit = 11
        ds.PixelRepresentation = 0  # unsigned
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.RescaleIntercept = str(int(data_min))
        ds.RescaleSlope = str(round((data_max - data_min) / 4095, 6))
        ds.WindowCenter = str(int((data_max + data_min) / 2))
        ds.WindowWidth = str(int(data_max - data_min))
        ds.PixelData = slice_data.astype(np.uint16).tobytes()

        filepath = series_dir / f"slice_{i:04d}.dcm"
        ds.save_as(str(filepath), write_like_original=False)

    print(f"  Written {len(slice_indices)} files to {series_dir}")
    return series_dir


def upload_to_orthanc(dicom_dir: Path, orthanc_url: str, username: str = "admin", password: str = "orthanc"):
    """Upload DICOM files to Orthanc via REST API.

    Args:
        dicom_dir: Directory containing .dcm files.
        orthanc_url: Orthanc REST API URL.
        username: Orthanc username.
        password: Orthanc password.
    """
    import requests

    dcm_files = sorted(dicom_dir.glob("*.dcm"))
    print(f"Uploading {len(dcm_files)} DICOM files to {orthanc_url}...")

    success = 0
    for dcm_path in dcm_files:
        with open(dcm_path, "rb") as f:
            resp = requests.post(
                f"{orthanc_url}/instances",
                data=f.read(),
                auth=(username, password),
                headers={"Content-Type": "application/dicom"},
                timeout=30,
            )
        if resp.status_code == 200:
            success += 1
        else:
            print(f"  Failed to upload {dcm_path.name}: {resp.status_code} {resp.text[:100]}")

    print(f"  Uploaded {success}/{len(dcm_files)} files successfully")


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic DICOM series from NIfTI samples")
    parser.add_argument("--upload", action="store_true", help="Upload generated DICOMs to Orthanc")
    parser.add_argument("--orthanc-url", default="http://localhost:8042", help="Orthanc REST API URL")
    parser.add_argument("--orthanc-user", default="admin", help="Orthanc username")
    parser.add_argument("--orthanc-pass", default="orthanc", help="Orthanc password")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: data/sample_dicoms/)")
    parser.add_argument("--max-slices", type=int, default=50, help="Max slices per series")
    args = parser.parse_args()

    sample_dir = PROJECT_ROOT / "data" / "sample_images"
    output_dir = Path(args.output_dir) if args.output_dir else PROJECT_ROOT / "data" / "sample_dicoms"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define samples to convert
    samples = [
        {
            "nifti": sample_dir / "sample_ct_chest.nii.gz",
            "patient_name": "DEMO^ChestCT",
            "patient_id": "DEMO-CHEST-001",
            "study_description": "CT Chest Demo",
            "series_description": "Axial Chest CT",
            "modality": "CT",
        },
        {
            "nifti": sample_dir / "sample_ct_head.nii.gz",
            "patient_name": "DEMO^HeadCT",
            "patient_id": "DEMO-HEAD-001",
            "study_description": "CT Head Demo",
            "series_description": "Axial Head CT",
            "modality": "CT",
        },
        {
            "nifti": sample_dir / "sample_brain_flair.nii.gz",
            "patient_name": "DEMO^BrainMRI",
            "patient_id": "DEMO-BRAIN-001",
            "study_description": "MRI Brain FLAIR Demo",
            "series_description": "Axial FLAIR",
            "modality": "MR",
        },
    ]

    generated_dirs = []
    for sample in samples:
        nifti_path = sample["nifti"]
        if not nifti_path.exists():
            print(f"Skipping {nifti_path.name} (not found)")
            continue
        series_dir = generate_dicom_series(
            nifti_path=nifti_path,
            output_dir=output_dir,
            patient_name=sample["patient_name"],
            patient_id=sample["patient_id"],
            study_description=sample["study_description"],
            series_description=sample["series_description"],
            modality=sample["modality"],
            max_slices=args.max_slices,
        )
        generated_dirs.append(series_dir)

    if args.upload and generated_dirs:
        print(f"\nUploading to Orthanc at {args.orthanc_url}...")
        for d in generated_dirs:
            upload_to_orthanc(d, args.orthanc_url, args.orthanc_user, args.orthanc_pass)

    print("\nDone!")


if __name__ == "__main__":
    main()
