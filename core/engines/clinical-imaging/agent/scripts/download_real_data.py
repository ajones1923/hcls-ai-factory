#!/usr/bin/env python3
"""Download real publicly available medical imaging data for validation.

Uses MedMNIST (CC BY 4.0) which provides real medical images from published
clinical datasets at 28x28 resolution (~10-50 MB per dataset). Images are
upscaled to 224x224 on save to match model input requirements.

Also generates full-resolution (1024x1024) synthetic chest X-ray images
with clinically realistic anatomy (lung fields, cardiac silhouette, rib
shadows, mediastinal structures) for validating the CXR workflow at native
clinical resolution without requiring multi-GB dataset downloads.

Datasets downloaded:
    - ChestMNIST: NIH ChestX-ray14 (14 pathology labels, multi-label)
    - PneumoniaMNIST: Guangzhou Women & Children's Medical Center (binary)
    - BreastMNIST: Breast ultrasound dataset (3-class)
    - OrganAMNIST: Abdominal CT organ classification (11-class)
    - Full-resolution CXR: Synthetic 1024x1024 PA chest X-rays

Also extracts a sample DICOM from pydicom's built-in test data.

Total download: ~50-100 MB (MedMNIST 28x28 default splits).
Images saved at 224x224 in: data/sample_images/
Full-resolution images saved at 1024x1024 in: data/sample_images/fullres/
Metadata stored in: data/sample_images/metadata.json

Author: Adam Jones
Date: 2026-02-27
License: Apache 2.0 (code), CC BY 4.0 (MedMNIST data)
"""

import json
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
SAMPLE_DIR = BASE_DIR / "data" / "sample_images"
FULLRES_DIR = SAMPLE_DIR / "fullres"
CACHE_DIR = SAMPLE_DIR / ".medmnist_cache"
SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
FULLRES_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Number of sample images to save per dataset
NUM_CHEST = 10
NUM_PNEUMONIA = 10
NUM_BREAST = 5
NUM_ORGAN = 5

# Output image resolution (matches CXR model input)
OUTPUT_SIZE = 224

# MedMNIST downloads at 28x28 (default, small file size)
# We upscale to 224x224 using bicubic interpolation on save
DOWNLOAD_SIZE = 28

# MedMNIST ChestMNIST label names (14 pathology labels, multi-label)
CHEST_LABEL_NAMES = [
    "Atelectasis",
    "Cardiomegaly",
    "Effusion",
    "Infiltration",
    "Mass",
    "Nodule",
    "Pneumonia",
    "Pneumothorax",
    "Consolidation",
    "Edema",
    "Emphysema",
    "Fibrosis",
    "Pleural Thickening",
    "Hernia",
]

# PneumoniaMNIST label names (binary)
PNEUMONIA_LABEL_NAMES = ["Normal", "Pneumonia"]

# BreastMNIST label names (3-class)
BREAST_LABEL_NAMES = ["Normal", "Benign", "Malignant"]

# OrganAMNIST label names (11-class)
ORGAN_LABEL_NAMES = [
    "Bladder",
    "Femur-L",
    "Femur-R",
    "Heart",
    "Kidney-L",
    "Kidney-R",
    "Liver",
    "Lung-L",
    "Lung-R",
    "Spleen",
    "Stomach",
]


def save_upscaled(img, filepath: Path) -> None:
    """Save a PIL image upscaled to OUTPUT_SIZE x OUTPUT_SIZE."""
    if not isinstance(img, Image.Image):
        img = Image.fromarray(np.array(img))
    img = img.resize((OUTPUT_SIZE, OUTPUT_SIZE), Image.BICUBIC)
    img.save(str(filepath))


def download_chest_xray() -> dict:
    """Download ChestMNIST (NIH ChestX-ray14) samples."""
    from medmnist import ChestMNIST

    print(f"\n{'='*60}")
    print("Downloading ChestMNIST (NIH ChestX-ray14)")
    print(f"{'='*60}")

    dataset = ChestMNIST(
        split="test", download=True, size=DOWNLOAD_SIZE,
        root=str(CACHE_DIR),
    )
    print(f"  Dataset size: {len(dataset)} test images")
    print(f"  Download size: {DOWNLOAD_SIZE}x{DOWNLOAD_SIZE}, output: {OUTPUT_SIZE}x{OUTPUT_SIZE}")
    print(f"  Labels: {len(CHEST_LABEL_NAMES)} pathologies (multi-label)")

    files = []

    # Select diverse samples: try to find images with different pathologies
    label_indices = {}
    for i in range(min(500, len(dataset))):
        _, label = dataset[i]
        label_array = np.array(label).flatten()
        for j, val in enumerate(label_array):
            if val == 1 and j not in label_indices:
                label_indices[j] = i

    # Build sample indices: prioritize diversity, then fill remaining
    selected = set()
    for label_idx in sorted(label_indices.keys()):
        if len(selected) >= NUM_CHEST:
            break
        selected.add(label_indices[label_idx])

    # Fill remaining with sequential indices
    idx = 0
    while len(selected) < NUM_CHEST and idx < len(dataset):
        selected.add(idx)
        idx += 1

    selected = sorted(selected)[:NUM_CHEST]

    for seq, i in enumerate(selected):
        img, label = dataset[i]
        label_array = np.array(label).flatten().tolist()

        filename = f"chest_xray_{seq:03d}.png"
        filepath = SAMPLE_DIR / filename
        save_upscaled(img, filepath)

        active_labels = [
            CHEST_LABEL_NAMES[j]
            for j in range(len(label_array))
            if j < len(CHEST_LABEL_NAMES) and label_array[j] == 1
        ]

        files.append({
            "filename": filename,
            "dataset_index": int(i),
            "label_vector": [int(x) for x in label_array],
            "label_names": active_labels if active_labels else ["No Finding"],
        })

        print(f"  Saved {filename} -> {active_labels if active_labels else ['No Finding']}")

    return {
        "source": "NIH ChestX-ray14 via MedMNIST v3",
        "license": "CC BY 4.0",
        "citation": "Wang et al., ChestX-ray8 (CVPR 2017); Yang et al., MedMNIST v2 (Sci Data 2023)",
        "num_samples": len(files),
        "total_dataset_size": len(dataset),
        "resolution": f"{OUTPUT_SIZE}x{OUTPUT_SIZE} (upscaled from {DOWNLOAD_SIZE}x{DOWNLOAD_SIZE})",
        "task": "multi-label classification (14 pathologies)",
        "label_names": CHEST_LABEL_NAMES,
        "files": files,
    }


def download_pneumonia() -> dict:
    """Download PneumoniaMNIST (pediatric CXR) samples."""
    from medmnist import PneumoniaMNIST

    print(f"\n{'='*60}")
    print("Downloading PneumoniaMNIST (Pediatric CXR)")
    print(f"{'='*60}")

    dataset = PneumoniaMNIST(
        split="test", download=True, size=DOWNLOAD_SIZE,
        root=str(CACHE_DIR),
    )
    print(f"  Dataset size: {len(dataset)} test images")

    files = []

    # Get balanced samples: some normal, some pneumonia
    normal_indices = []
    pneumonia_indices = []

    for i in range(min(300, len(dataset))):
        _, label = dataset[i]
        label_val = int(np.array(label).flatten()[0])
        if label_val == 0 and len(normal_indices) < NUM_PNEUMONIA // 2 + 1:
            normal_indices.append(i)
        elif label_val == 1 and len(pneumonia_indices) < NUM_PNEUMONIA // 2 + 1:
            pneumonia_indices.append(i)
        if len(normal_indices) + len(pneumonia_indices) >= NUM_PNEUMONIA:
            break

    selected = sorted(normal_indices + pneumonia_indices)[:NUM_PNEUMONIA]

    for seq, i in enumerate(selected):
        img, label = dataset[i]
        label_val = int(np.array(label).flatten()[0])

        filename = f"pneumonia_{seq:03d}.png"
        filepath = SAMPLE_DIR / filename
        save_upscaled(img, filepath)

        label_name = PNEUMONIA_LABEL_NAMES[label_val]
        files.append({
            "filename": filename,
            "dataset_index": int(i),
            "label": label_val,
            "label_name": label_name,
        })

        print(f"  Saved {filename} -> {label_name}")

    return {
        "source": "Guangzhou Women and Children's Medical Center via MedMNIST v3",
        "license": "CC BY 4.0",
        "citation": "Kermany et al., Cell 2018; Yang et al., MedMNIST v2 (Sci Data 2023)",
        "num_samples": len(files),
        "total_dataset_size": len(dataset),
        "resolution": f"{OUTPUT_SIZE}x{OUTPUT_SIZE} (upscaled from {DOWNLOAD_SIZE}x{DOWNLOAD_SIZE})",
        "task": "binary classification (normal vs pneumonia)",
        "label_names": PNEUMONIA_LABEL_NAMES,
        "files": files,
    }


def download_breast() -> dict:
    """Download BreastMNIST (breast ultrasound) samples."""
    from medmnist import BreastMNIST

    print(f"\n{'='*60}")
    print("Downloading BreastMNIST (Breast Ultrasound)")
    print(f"{'='*60}")

    dataset = BreastMNIST(
        split="test", download=True, size=DOWNLOAD_SIZE,
        root=str(CACHE_DIR),
    )
    print(f"  Dataset size: {len(dataset)} test images")

    files = []

    # Try to get representative samples from each class
    class_indices = {0: [], 1: [], 2: []}
    for i in range(len(dataset)):
        _, label = dataset[i]
        label_val = int(np.array(label).flatten()[0])
        if label_val in class_indices and len(class_indices[label_val]) < 2:
            class_indices[label_val].append(i)
        if sum(len(v) for v in class_indices.values()) >= NUM_BREAST:
            break

    # If we did not find enough classes, just take sequential samples
    selected = sorted(
        [idx for indices in class_indices.values() for idx in indices]
    )
    if len(selected) < NUM_BREAST:
        for i in range(len(dataset)):
            if i not in selected:
                selected.append(i)
            if len(selected) >= NUM_BREAST:
                break
        selected = sorted(selected)
    selected = selected[:NUM_BREAST]

    for seq, i in enumerate(selected):
        img, label = dataset[i]
        label_val = int(np.array(label).flatten()[0])

        filename = f"breast_us_{seq:03d}.png"
        filepath = SAMPLE_DIR / filename
        save_upscaled(img, filepath)

        label_name = BREAST_LABEL_NAMES[label_val] if label_val < len(BREAST_LABEL_NAMES) else f"Class_{label_val}"
        files.append({
            "filename": filename,
            "dataset_index": int(i),
            "label": label_val,
            "label_name": label_name,
        })

        print(f"  Saved {filename} -> {label_name}")

    return {
        "source": "Breast Ultrasound Images Dataset (BUSI) via MedMNIST v3",
        "license": "CC BY 4.0",
        "citation": "Al-Dhabyani et al., Data in Brief 2020; Yang et al., MedMNIST v2 (Sci Data 2023)",
        "num_samples": len(files),
        "total_dataset_size": len(dataset),
        "resolution": f"{OUTPUT_SIZE}x{OUTPUT_SIZE} (upscaled from {DOWNLOAD_SIZE}x{DOWNLOAD_SIZE})",
        "task": "3-class classification (normal, benign, malignant)",
        "label_names": BREAST_LABEL_NAMES,
        "files": files,
    }


def download_organ() -> dict:
    """Download OrganAMNIST (abdominal CT axial slices) samples."""
    from medmnist import OrganAMNIST

    print(f"\n{'='*60}")
    print("Downloading OrganAMNIST (Abdominal CT)")
    print(f"{'='*60}")

    dataset = OrganAMNIST(
        split="test", download=True, size=DOWNLOAD_SIZE,
        root=str(CACHE_DIR),
    )
    print(f"  Dataset size: {len(dataset)} test images")

    files = []

    # Try to get diverse organ classes
    class_indices = {}
    for i in range(min(500, len(dataset))):
        _, label = dataset[i]
        label_val = int(np.array(label).flatten()[0])
        if label_val not in class_indices:
            class_indices[label_val] = i
        if len(class_indices) >= NUM_ORGAN:
            break

    selected = sorted(class_indices.values())[:NUM_ORGAN]

    for seq, i in enumerate(selected):
        img, label = dataset[i]
        label_val = int(np.array(label).flatten()[0])

        filename = f"organ_ct_{seq:03d}.png"
        filepath = SAMPLE_DIR / filename
        save_upscaled(img, filepath)

        label_name = ORGAN_LABEL_NAMES[label_val] if label_val < len(ORGAN_LABEL_NAMES) else f"Organ_{label_val}"
        files.append({
            "filename": filename,
            "dataset_index": int(i),
            "label": label_val,
            "label_name": label_name,
        })

        print(f"  Saved {filename} -> {label_name}")

    return {
        "source": "Liver Tumor Segmentation Benchmark (LiTS) via MedMNIST v3",
        "license": "CC BY 4.0",
        "citation": "Bilic et al., Medical Image Analysis 2023; Yang et al., MedMNIST v2 (Sci Data 2023)",
        "num_samples": len(files),
        "total_dataset_size": len(dataset),
        "resolution": f"{OUTPUT_SIZE}x{OUTPUT_SIZE} (upscaled from {DOWNLOAD_SIZE}x{DOWNLOAD_SIZE})",
        "task": "11-class classification (organ identification from axial CT slices)",
        "label_names": ORGAN_LABEL_NAMES,
        "files": files,
    }


def download_dicom_sample() -> dict:
    """Extract a sample DICOM from pydicom's built-in test files."""
    print(f"\n{'='*60}")
    print("Extracting DICOM sample from pydicom test data")
    print(f"{'='*60}")

    try:
        import pydicom
        from pydicom.data import get_testdata_file

        # pydicom ships with several test DICOM files
        test_file = get_testdata_file("CT_small.dcm")
        ds = pydicom.dcmread(test_file)

        filename = "dicom_ct_sample.dcm"
        filepath = SAMPLE_DIR / filename
        ds.save_as(str(filepath))

        # Also save as PNG for quick viewing
        pixel_array = ds.pixel_array.astype(np.float32)
        if pixel_array.max() > 0:
            normalized = ((pixel_array - pixel_array.min()) /
                          (pixel_array.max() - pixel_array.min()) * 255).astype(np.uint8)
        else:
            normalized = np.zeros_like(pixel_array, dtype=np.uint8)

        png_filename = "dicom_ct_sample_preview.png"
        png_filepath = SAMPLE_DIR / png_filename
        pil_img = Image.fromarray(normalized)
        pil_img = pil_img.resize((OUTPUT_SIZE, OUTPUT_SIZE), Image.BICUBIC)
        pil_img.save(str(png_filepath))

        info = {
            "filename": filename,
            "preview_png": png_filename,
            "modality": str(getattr(ds, "Modality", "Unknown")),
            "patient_id": str(getattr(ds, "PatientID", "Anonymous")),
            "rows": int(getattr(ds, "Rows", 0)),
            "columns": int(getattr(ds, "Columns", 0)),
            "bits_stored": int(getattr(ds, "BitsStored", 0)),
            "pixel_spacing": [float(x) for x in ds.PixelSpacing] if hasattr(ds, "PixelSpacing") else [],
            "slice_thickness": float(ds.SliceThickness) if hasattr(ds, "SliceThickness") else None,
        }

        print(f"  Saved {filename} ({info['modality']}, {info['rows']}x{info['columns']})")
        print(f"  Saved {png_filename} (preview)")

        return {
            "source": "pydicom test data (public domain CT slice)",
            "license": "MIT (pydicom), public domain (image data)",
            "files": [info],
        }

    except Exception as e:
        print(f"  Warning: Could not extract DICOM sample: {e}")
        return {"source": "pydicom test data", "error": str(e), "files": []}


# ---------------------------------------------------------------------------
# Full-Resolution Synthetic CXR Generation
# ---------------------------------------------------------------------------

# Publicly accessible CXR sample URLs to try before falling back to
# synthetic generation.  These are small individual files that do not
# require authentication.
_PUBLIC_CXR_URLS = [
    # NIH Clinical Center sample CXR images hosted on Box.com (public links)
    "https://nihcc.box.com/shared/static/vfk49d74nhbxq3nqjg0900w5nvkorp5c.png",
    "https://nihcc.box.com/shared/static/i28rlmbvmfjbl8p2n3ril0pez8ffoosh.png",
    "https://nihcc.box.com/shared/static/f1t00wrtdk94satdfb9olcolqx20z2jp.png",
    # OpenI / Indiana University CXR images (direct PNG links)
    "https://openi.nlm.nih.gov/imgs/collections/NLMCXR_png/CXR1_1_IM-0001-3001.png",
    "https://openi.nlm.nih.gov/imgs/collections/NLMCXR_png/CXR1_1_IM-0002-1001.png",
]


def _try_download_public_cxr(output_dir: Path, count: int) -> list:
    """Attempt to download full-resolution CXR images from public URLs.

    Returns list of dicts with metadata for successfully downloaded images.
    """
    downloaded = []
    for idx, url in enumerate(_PUBLIC_CXR_URLS):
        if len(downloaded) >= count:
            break
        filename = f"fullres_cxr_real_{idx:03d}.png"
        filepath = output_dir / filename
        if filepath.exists():
            print(f"  Already exists: {filename}")
            img = Image.open(str(filepath))
            downloaded.append({
                "filename": filename,
                "source_url": url,
                "source": "public_download",
                "resolution": f"{img.size[0]}x{img.size[1]}",
                "width": img.size[0],
                "height": img.size[1],
            })
            continue

        try:
            print(f"  Trying URL: {url[:80]}...")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = resp.read()
            # Validate that it is actually an image
            import io
            img = Image.open(io.BytesIO(data))
            if img.size[0] < 256 or img.size[1] < 256:
                print(f"    Skipped: too small ({img.size[0]}x{img.size[1]})")
                continue
            img.save(str(filepath))
            print(f"    Downloaded {filename} ({img.size[0]}x{img.size[1]})")
            downloaded.append({
                "filename": filename,
                "source_url": url,
                "source": "public_download",
                "resolution": f"{img.size[0]}x{img.size[1]}",
                "width": img.size[0],
                "height": img.size[1],
            })
        except Exception as e:
            print(f"    Failed: {e}")
            continue

    return downloaded


def _generate_synthetic_cxr(
    size: int = 1024,
    seed: int = 0,
    pathology: str = "normal",
) -> np.ndarray:
    """Generate a clinically realistic synthetic PA chest X-ray image.

    Creates an image with the following anatomical structures:
        - Background soft-tissue density (uniform gray)
        - Lung fields (dark ovoid regions with slight texture)
        - Cardiac silhouette (bright region in left-lower hemithorax)
        - Mediastinal structures (bright central vertical band)
        - Rib shadows (periodic horizontal bright bands with curvature)
        - Diaphragm domes (bright curved inferior border)
        - Spine shadow (faint midline vertical density)
        - Optional pathology overlays: consolidation, effusion, cardiomegaly

    Args:
        size: Output image dimension (square).
        seed: Random seed for reproducibility.
        pathology: One of "normal", "consolidation", "effusion",
            "cardiomegaly", "pneumothorax".

    Returns:
        uint8 numpy array of shape (size, size), values 0-255.
    """
    rng = np.random.RandomState(seed)
    H = W = size

    # Coordinate grids normalized to [-1, 1]
    y_grid, x_grid = np.mgrid[0:H, 0:W].astype(np.float64)
    y_norm = (y_grid / H) * 2 - 1  # -1 (top) to +1 (bottom)
    x_norm = (x_grid / W) * 2 - 1  # -1 (left) to +1 (right)

    # Start with a uniform soft-tissue background (moderate gray)
    canvas = np.full((H, W), 160.0, dtype=np.float64)

    # ── Mediastinum (bright central vertical band) ──
    mediastinum_width = 0.15
    mediastinum_mask = np.exp(-(x_norm ** 2) / (2 * mediastinum_width ** 2))
    canvas += mediastinum_mask * 40

    # ── Spine shadow (faint midline) ──
    spine_width = 0.04
    spine_mask = np.exp(-(x_norm ** 2) / (2 * spine_width ** 2))
    canvas += spine_mask * 15

    # ── Lung fields (dark ovoid regions) ──
    # Right lung (patient's right = image left in PA view)
    right_lung_cx, right_lung_cy = -0.30, -0.05
    right_lung_sx, right_lung_sy = 0.25, 0.40
    right_lung = np.exp(
        -((x_norm - right_lung_cx) ** 2) / (2 * right_lung_sx ** 2)
        - ((y_norm - right_lung_cy) ** 2) / (2 * right_lung_sy ** 2)
    )

    # Left lung (patient's left = image right in PA view)
    left_lung_cx, left_lung_cy = 0.30, -0.05
    left_lung_sx, left_lung_sy = 0.23, 0.38
    left_lung = np.exp(
        -((x_norm - left_lung_cx) ** 2) / (2 * left_lung_sx ** 2)
        - ((y_norm - left_lung_cy) ** 2) / (2 * left_lung_sy ** 2)
    )

    lung_mask = np.maximum(right_lung, left_lung)
    canvas -= lung_mask * 90  # Lungs are darker (more air = more X-ray penetration)

    # ── Lung texture (subtle noise within lung fields) ──
    lung_noise = rng.randn(H, W) * 5
    lung_noise = ndimage.gaussian_filter(lung_noise, sigma=size / 80)
    canvas += lung_noise * lung_mask

    # ── Cardiac silhouette (bright mass, left-lower) ──
    heart_cx, heart_cy = 0.05, 0.15  # Slightly left of center, lower
    heart_sx, heart_sy = 0.18, 0.20
    if pathology == "cardiomegaly":
        heart_sx *= 1.5
        heart_sy *= 1.3
    heart_mask = np.exp(
        -((x_norm - heart_cx) ** 2) / (2 * heart_sx ** 2)
        - ((y_norm - heart_cy) ** 2) / (2 * heart_sy ** 2)
    )
    canvas += heart_mask * 70

    # ── Diaphragm domes (bright curved band at lung bases) ──
    # Right dome (slightly higher than left)
    for side, cx, dome_y, amplitude in [(-0.30, -0.30, 0.42, 0.08), (0.30, 0.30, 0.45, 0.07)]:
        diaphragm_curve = dome_y + amplitude * np.cos(np.pi * (x_norm - cx) / 0.5)
        diaphragm_band = np.exp(-((y_norm - diaphragm_curve) ** 2) / (2 * 0.03 ** 2))
        # Only apply where lungs are (mask with x-position)
        x_mask = np.exp(-((x_norm - cx) ** 2) / (2 * 0.25 ** 2))
        canvas += diaphragm_band * x_mask * 50

    # ── Rib shadows (periodic curved bright bands) ──
    n_ribs = 10
    for rib_idx in range(n_ribs):
        rib_y_center = -0.60 + rib_idx * 0.10  # Spaced vertically
        # Ribs curve downward laterally
        rib_curvature = 0.03 * (rib_idx + 1)
        rib_curve = rib_y_center + rib_curvature * (x_norm ** 2)
        rib_band = np.exp(-((y_norm - rib_curve) ** 2) / (2 * 0.008 ** 2))
        # Ribs are wider laterally, narrow at spine
        lateral_weight = np.clip(np.abs(x_norm) * 2, 0.2, 1.0)
        # Only in lung field region
        rib_intensity = 12 + rng.rand() * 5
        canvas += rib_band * lateral_weight * rib_intensity * lung_mask

    # ── Clavicles (bright diagonal bands at top) ──
    for sign in [-1, 1]:
        clav_slope = sign * 0.15
        clav_y = -0.55 + clav_slope * x_norm
        clav_band = np.exp(-((y_norm - clav_y) ** 2) / (2 * 0.012 ** 2))
        clav_x_mask = np.exp(-((x_norm - sign * 0.25) ** 2) / (2 * 0.20 ** 2))
        canvas += clav_band * clav_x_mask * 25

    # ── Scapulae (faint lateral bright regions) ──
    for sign in [-1, 1]:
        scap_cx = sign * 0.55
        scap_cy = -0.15
        scap_mask = np.exp(
            -((x_norm - scap_cx) ** 2) / (2 * 0.10 ** 2)
            - ((y_norm - scap_cy) ** 2) / (2 * 0.25 ** 2)
        )
        canvas += scap_mask * 20

    # ── Pathology overlays ──
    if pathology == "consolidation":
        # Dense opacity in right lower lobe
        consol_cx, consol_cy = -0.25, 0.20
        consol_sx, consol_sy = 0.12, 0.10
        consol_mask = np.exp(
            -((x_norm - consol_cx) ** 2) / (2 * consol_sx ** 2)
            - ((y_norm - consol_cy) ** 2) / (2 * consol_sy ** 2)
        )
        canvas += consol_mask * 60
        # Add air bronchograms (small linear lucencies within consolidation)
        for _ in range(3):
            ab_angle = rng.uniform(-0.3, 0.3)
            ab_y = consol_cy + rng.uniform(-0.05, 0.05)
            ab_line = np.exp(-((y_norm - ab_y - ab_angle * x_norm) ** 2) / (2 * 0.004 ** 2))
            ab_x_mask = np.exp(-((x_norm - consol_cx) ** 2) / (2 * 0.08 ** 2))
            canvas -= ab_line * ab_x_mask * consol_mask * 25

    elif pathology == "effusion":
        # Meniscus sign: fluid layering at right costophrenic angle
        effusion_level = 0.35  # Upper boundary of effusion
        effusion_mask = np.clip((y_norm - effusion_level) / 0.15, 0, 1)
        # Only in right hemithorax
        right_mask = np.exp(-((x_norm + 0.30) ** 2) / (2 * 0.22 ** 2))
        # Meniscus curve (fluid rises at periphery)
        meniscus = 0.05 * np.exp(-((x_norm + 0.30) ** 2) / (2 * 0.10 ** 2))
        effusion_curved = np.clip((y_norm - effusion_level + meniscus) / 0.15, 0, 1)
        canvas += effusion_curved * right_mask * 55

    elif pathology == "pneumothorax":
        # Visceral pleural line with absent lung markings laterally
        ptx_side = -1  # Right-sided
        pleural_line_x = ptx_side * 0.40
        # Thin bright pleural line
        pl_mask = np.exp(-((x_norm - pleural_line_x) ** 2) / (2 * 0.005 ** 2))
        pl_y_mask = np.clip(1 - (y_norm + 0.4) / 0.6, 0, 1)  # Upper lung
        canvas += pl_mask * pl_y_mask * 30
        # Absent lung markings lateral to pleural line
        lateral_mask = np.clip((ptx_side * x_norm - ptx_side * pleural_line_x) / 0.15, 0, 1)
        canvas += lateral_mask * pl_y_mask * 20  # Hyperlucent

    # ── Global noise and finishing ──
    global_noise = rng.randn(H, W) * 3
    global_noise = ndimage.gaussian_filter(global_noise, sigma=size / 200)
    canvas += global_noise

    # Vignette (slight darkening at edges, simulating collimation)
    vignette = 1.0 - 0.15 * (x_norm ** 2 + y_norm ** 2)
    canvas *= vignette

    # Clamp and convert to uint8
    canvas = np.clip(canvas, 0, 255).astype(np.uint8)

    return canvas


def download_fullres_cxr(output_dir: Path = FULLRES_DIR, count: int = 5) -> dict:
    """Download or generate full-resolution (1024x1024) chest X-ray images.

    Strategy (in order of preference):
        1. Try to download real CXR images from publicly accessible URLs
           (NIH Clinical Center, OpenI/Indiana University).
        2. Fall back to generating high-quality synthetic 1024x1024 CXR
           images with realistic anatomical structures.

    All images are saved as 8-bit grayscale PNGs in output_dir.

    Args:
        output_dir: Directory to save images into.
        count: Number of images to produce.

    Returns:
        Metadata dict compatible with the existing metadata.json schema.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print("Full-Resolution CXR Images (1024x1024)")
    print(f"{'='*60}")
    print(f"  Output directory: {output_dir}")
    print(f"  Target count: {count}")

    # ── Step 1: Try public downloads ──
    print("\n  Step 1: Attempting public CXR downloads...")
    downloaded = _try_download_public_cxr(output_dir, count)
    print(f"  Downloaded {len(downloaded)} real CXR images")

    # ── Step 2: Generate synthetic CXRs for any remaining ──
    n_remaining = count - len(downloaded)
    generated = []

    if n_remaining > 0:
        print(f"\n  Step 2: Generating {n_remaining} synthetic 1024x1024 CXR images...")

        # Create a variety of pathologies
        pathologies = ["normal", "consolidation", "effusion", "cardiomegaly", "pneumothorax"]

        for i in range(n_remaining):
            pathology = pathologies[i % len(pathologies)]
            seed = 42 + i

            print(f"    Generating fullres_cxr_synth_{i:03d}.png (pathology={pathology})...")
            img_array = _generate_synthetic_cxr(size=1024, seed=seed, pathology=pathology)

            filename = f"fullres_cxr_synth_{i:03d}.png"
            filepath = output_dir / filename
            img = Image.fromarray(img_array, mode="L")
            img.save(str(filepath))

            generated.append({
                "filename": filename,
                "source": "synthetic",
                "pathology": pathology,
                "seed": seed,
                "resolution": "1024x1024",
                "width": 1024,
                "height": 1024,
            })

            print(f"      Saved {filename} ({img_array.shape[0]}x{img_array.shape[1]}, "
                  f"range=[{img_array.min()}, {img_array.max()}])")

    all_files = downloaded + generated

    # Disk usage
    total_bytes = sum(
        (output_dir / f["filename"]).stat().st_size
        for f in all_files
        if (output_dir / f["filename"]).exists()
    )

    print(f"\n  Summary:")
    print(f"    Real downloads: {len(downloaded)}")
    print(f"    Synthetic generated: {len(generated)}")
    print(f"    Total images: {len(all_files)}")
    print(f"    Total disk usage: {total_bytes / 1024:.1f} KB")

    return {
        "source": "Synthetic full-resolution CXR (1024x1024) with realistic anatomy",
        "license": "Apache 2.0 (synthetic), CC BY 4.0 (any downloaded real images)",
        "citation": "Synthetic CXR generated with Gaussian anatomical models",
        "num_samples": len(all_files),
        "resolution": "1024x1024 (native full-resolution)",
        "task": "CXR multi-label classification at clinical resolution",
        "generation_details": (
            "Anatomical structures: lung fields, cardiac silhouette, ribs, "
            "mediastinum, diaphragm, clavicles, scapulae. "
            "Pathology overlays: consolidation, effusion, cardiomegaly, pneumothorax."
        ),
        "files": all_files,
        "n_real_downloads": len(downloaded),
        "n_synthetic": len(generated),
    }


def main():
    """Download all datasets and write metadata JSON."""
    print("=" * 60)
    print("MedMNIST Real Medical Imaging Data Download")
    print("=" * 60)
    print(f"Output directory: {SAMPLE_DIR}")
    print(f"Download resolution: {DOWNLOAD_SIZE}x{DOWNLOAD_SIZE}")
    print(f"Output resolution: {OUTPUT_SIZE}x{OUTPUT_SIZE}")
    print()

    metadata = {}

    # Download each dataset
    metadata["chest_xray"] = download_chest_xray()
    metadata["pneumonia"] = download_pneumonia()
    metadata["breast_ultrasound"] = download_breast()
    metadata["organ_ct"] = download_organ()
    metadata["dicom_sample"] = download_dicom_sample()
    metadata["fullres_cxr"] = download_fullres_cxr()

    # Write metadata
    metadata_path = SAMPLE_DIR / "metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"\nMetadata written to: {metadata_path}")

    # Summary
    total_images = sum(
        d.get("num_samples", len(d.get("files", [])))
        for d in metadata.values()
    )
    print(f"\n{'='*60}")
    print("DOWNLOAD SUMMARY")
    print(f"{'='*60}")
    print(f"  Total images saved: {total_images}")
    print(f"  Datasets: {len(metadata)}")
    for name, data in metadata.items():
        n = data.get("num_samples", len(data.get("files", [])))
        src = data.get("source", "unknown")
        print(f"    {name}: {n} images ({src})")
    print(f"\n  All files in: {SAMPLE_DIR}")
    print(f"  Metadata: {metadata_path}")

    # Disk usage
    total_bytes = 0
    for pattern in ("*.png", "*.dcm"):
        for f in SAMPLE_DIR.glob(pattern):
            total_bytes += f.stat().st_size
    print(f"  Total image size on disk: {total_bytes / 1024 / 1024:.1f} MB")

    # Cache size
    cache_bytes = sum(
        f.stat().st_size
        for f in CACHE_DIR.rglob("*")
        if f.is_file()
    )
    print(f"  MedMNIST cache size: {cache_bytes / 1024 / 1024:.1f} MB")

    return 0


if __name__ == "__main__":
    sys.exit(main())
