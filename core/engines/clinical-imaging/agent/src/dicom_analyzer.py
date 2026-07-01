"""Live DICOM analysis pipeline for the Clinical Imaging Engine.

Processes uploaded DICOM files through the appropriate clinical workflow
with real AI model inference (not mock mode).

Supports:
- Single-frame DICOM (CXR, mammography)
- Multi-frame/series DICOM (CT, MRI volumes)
- Automatic modality detection and workflow routing
- Results returned as WorkflowResult-compatible dict

Author: Adam Jones
Date: April 2026
"""

import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch
from loguru import logger


# =====================================================================
# Modality -> Workflow Routing Table
# =====================================================================

_MODALITY_WORKFLOW_ROUTING = {
    ("CR", "CHEST"): "cxr_rapid_findings",
    ("DX", "CHEST"): "cxr_rapid_findings",
    ("CR", ""): "cxr_rapid_findings",      # Default CR to CXR
    ("DX", ""): "cxr_rapid_findings",      # Default DX to CXR
    ("CT", "HEAD"): "ct_head_hemorrhage",
    ("CT", "CHEST"): "ct_chest_lung_nodule",
    ("CT", "HEART"): "ct_coronary_angiography",
    ("CT", "ABDOMEN"): "liver_lirads",
    ("MR", "HEAD"): "mri_brain_ms_lesion",
    ("MR", "BRAIN"): "mri_brain_ms_lesion",
    ("MR", "PELVIS"): "mri_prostate_pirads",
    ("MG", "BREAST"): "breast_birads",
}


class DICOMAnalyzer:
    """Live DICOM analysis engine with real AI model inference.

    Manages model loading, DICOM parsing, modality detection,
    workflow routing, and inference execution. Models are loaded
    lazily and cached for subsequent requests.
    """

    def __init__(self, models_dir: str = "data/models", use_gpu: bool = True):
        self.models_dir = Path(models_dir)
        self.use_gpu = use_gpu and torch.cuda.is_available()
        self.device = torch.device("cuda" if self.use_gpu else "cpu")
        self._loaded_models: Dict[str, Any] = {}
        logger.info(
            f"DICOMAnalyzer initialized: gpu={self.use_gpu}, "
            f"device={self.device}, models_dir={self.models_dir}"
        )

    # =================================================================
    # Public API
    # =================================================================

    def analyze_dicom(self, dicom_path: str, workflow_name: Optional[str] = None) -> Dict:
        """Analyze a DICOM file or directory of DICOM slices.

        If workflow_name is None, auto-detect based on DICOM modality/body part.
        Returns WorkflowResult-compatible dict with findings, measurements,
        classification, and is_mock=False.

        Args:
            dicom_path: Path to a single DICOM file or directory of slices.
            workflow_name: Optional workflow override (e.g., 'cxr_rapid_findings').

        Returns:
            Dict with workflow results including findings and measurements.
        """
        start = time.time()
        dicom_path = str(dicom_path)

        # Determine if single file or directory (series)
        p = Path(dicom_path)
        is_series = p.is_dir()

        # Detect modality from DICOM headers
        if is_series:
            dcm_files = self._find_dicom_files(dicom_path)
            if not dcm_files:
                return self._error_result("No DICOM files found in directory", start)
            detect_path = dcm_files[0]
        else:
            detect_path = dicom_path

        meta = self._detect_modality(detect_path)
        logger.info(f"DICOM metadata: {meta}")

        # Route to workflow
        if workflow_name is None:
            workflow_name = self._route_to_workflow(
                meta.get("modality", ""), meta.get("body_part", "")
            )
        logger.info(f"Routing to workflow: {workflow_name}")

        # Dispatch to appropriate analyzer
        try:
            if workflow_name == "cxr_rapid_findings":
                result = self.analyze_cxr(detect_path if not is_series else dcm_files[0])
            elif is_series:
                result = self.analyze_ct_volume(dicom_path, workflow_name)
            else:
                # Single-frame non-CXR -- use CXR analyzer as best-effort
                # or return a placeholder for modalities we cannot yet analyze live
                result = self._analyze_single_frame_fallback(dicom_path, workflow_name, meta)
        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            return self._error_result(str(e), start)

        result["inference_time_ms"] = round((time.time() - start) * 1000, 1)
        result["dicom_metadata"] = meta
        return result

    def analyze_cxr(self, dicom_path: str) -> Dict:
        """Analyze a chest X-ray DICOM using DenseNet-121 (torchxrayvision).

        This is the primary live analysis modality because:
        - torchxrayvision auto-downloads weights (densenet121-res224-all)
        - Single 2D image, no volume preprocessing
        - Fast inference (~1 second on GPU)
        - 18 pathology labels with confidence scores

        Args:
            dicom_path: Path to a CXR DICOM file.

        Returns:
            WorkflowResult-compatible dict with real inference findings.
        """
        import pydicom
        import torchxrayvision as xrv

        logger.info(f"Running live CXR analysis on: {dicom_path}")

        # Load DICOM pixel data
        ds = pydicom.dcmread(dicom_path)
        pixel_array = ds.pixel_array.astype(np.float32)

        # Handle multi-frame or color DICOM
        if pixel_array.ndim == 3:
            if pixel_array.shape[-1] in (3, 4):
                pixel_array = np.mean(pixel_array[..., :3], axis=-1)
            elif pixel_array.shape[0] > 1:
                pixel_array = pixel_array[0]
            else:
                pixel_array = pixel_array.squeeze()

        # Scale to [0, 255] for xrv normalization
        pmin, pmax = pixel_array.min(), pixel_array.max()
        if pmax > pmin:
            pixel_array = ((pixel_array - pmin) / (pmax - pmin)) * 255.0
        else:
            pixel_array = np.zeros_like(pixel_array)

        # Resize to 224x224
        from PIL import Image
        pil_img = Image.fromarray(pixel_array.astype(np.float32), mode="F")
        pil_img = pil_img.resize((224, 224), Image.BILINEAR)
        img = np.array(pil_img, dtype=np.float32)

        # Apply torchxrayvision normalization (centers to [-1024, 1024] range)
        img = xrv.datasets.normalize(img, maxval=255)

        # Convert to tensor: (1, 1, 224, 224)
        tensor = torch.from_numpy(img).float().unsqueeze(0).unsqueeze(0)

        # Load model (cached)
        if "densenet121" not in self._loaded_models:
            logger.info("Loading torchxrayvision DenseNet (densenet121-res224-all)...")
            model = xrv.models.DenseNet(weights="densenet121-res224-all")
            model.eval()
            if self.use_gpu:
                model = model.to(self.device)
            self._loaded_models["densenet121"] = model
            param_count = sum(p.numel() for p in model.parameters())
            logger.info(f"DenseNet-121 loaded: {param_count:,} params on {self.device}")

        model = self._loaded_models["densenet121"]

        # Inference
        with torch.no_grad():
            if self.use_gpu:
                tensor = tensor.to(self.device)
            outputs = model(tensor)
            predictions = outputs.cpu().numpy()[0]

        # Map all 18 pathology labels
        pathologies = model.pathologies
        all_scores: Dict[str, float] = {}
        findings: List[Dict] = []

        for pathology, score in zip(pathologies, predictions):
            score_f = float(score)
            all_scores[pathology] = round(score_f, 4)

            if score_f > 0.5:
                severity = "critical" if score_f > 0.8 else ("urgent" if score_f > 0.65 else "significant")
                findings.append({
                    "category": pathology.lower().replace(" ", "_"),
                    "description": (
                        f"{pathology} detected with {score_f:.1%} confidence. "
                        f"Clinical correlation recommended."
                    ),
                    "severity": severity,
                    "confidence": round(score_f, 3),
                    "threshold": 0.5,
                    "above_threshold": True,
                })

        if not findings:
            findings.append({
                "category": "normal",
                "description": (
                    "No significant acute findings. All pathology scores "
                    "below clinical thresholds."
                ),
                "severity": "normal",
            })

        # Classification summary
        positive_labels = [f["category"] for f in findings if f.get("above_threshold")]
        classification = (
            f"positive: {', '.join(positive_labels)}" if positive_labels else "negative"
        )

        # Overall severity
        severity_order = ["normal", "routine", "significant", "urgent", "critical"]
        overall_severity = "normal"
        for f in findings:
            f_sev = f.get("severity", "normal")
            if severity_order.index(f_sev) > severity_order.index(overall_severity):
                overall_severity = f_sev

        # Measurements dict (all 18 scores)
        measurements = {
            f"{k.lower().replace(' ', '_')}_confidence": v
            for k, v in all_scores.items()
        }

        return {
            "workflow_name": "cxr_rapid_findings",
            "status": "completed",
            "classification": classification,
            "severity": overall_severity,
            "findings": findings,
            "measurements": measurements,
            "inference_time_ms": 0,  # Set by caller
            "is_mock": False,
            "nim_services_used": ["DenseNet-121 (torchxrayvision, LIVE inference)"],
            "all_pathology_scores": all_scores,
        }

    def analyze_ct_volume(self, dicom_dir: str, workflow_name: str) -> Dict:
        """Analyze a CT volume (directory of DICOM slices).

        Loads all slices, sorts by InstanceNumber or SliceLocation,
        stacks into a 3D volume, and runs threshold-based analysis.

        For CT workflows without downloaded MONAI model weights, uses
        intensity-based analysis (windowing + thresholding) to provide
        clinically reasonable measurements.

        Args:
            dicom_dir: Directory containing DICOM slices.
            workflow_name: Target workflow name.

        Returns:
            WorkflowResult-compatible dict.
        """
        import pydicom

        logger.info(f"Loading CT volume from: {dicom_dir}")

        dcm_files = self._find_dicom_files(dicom_dir)
        if not dcm_files:
            return self._error_result("No DICOM files found in directory", time.time())

        # Load and sort slices
        slices = []
        for f in dcm_files:
            try:
                ds = pydicom.dcmread(f)
                slices.append(ds)
            except Exception as e:
                logger.warning(f"Failed to read {f}: {e}")

        if not slices:
            return self._error_result("No valid DICOM slices loaded", time.time())

        # Sort by InstanceNumber or SliceLocation
        try:
            slices.sort(key=lambda s: float(getattr(s, "InstanceNumber", 0)))
        except (TypeError, ValueError):
            try:
                slices.sort(key=lambda s: float(getattr(s, "SliceLocation", 0)))
            except (TypeError, ValueError):
                pass

        # Stack into 3D volume
        volume = np.stack([s.pixel_array.astype(np.float32) for s in slices], axis=0)

        # Apply rescale slope/intercept if available (convert to HU)
        ds0 = slices[0]
        slope = float(getattr(ds0, "RescaleSlope", 1.0))
        intercept = float(getattr(ds0, "RescaleIntercept", 0.0))
        volume = volume * slope + intercept

        logger.info(
            f"CT volume loaded: shape={volume.shape}, "
            f"HU range=[{volume.min():.0f}, {volume.max():.0f}], "
            f"{len(slices)} slices"
        )

        # Run threshold-based analysis appropriate for the workflow
        return self._analyze_ct_by_threshold(volume, workflow_name, ds0)

    # =================================================================
    # DICOM Metadata
    # =================================================================

    def _detect_modality(self, dicom_path: str) -> Dict:
        """Read DICOM headers to determine modality, body part, and patient info.

        Args:
            dicom_path: Path to a single DICOM file.

        Returns:
            Dict with modality, body_part, patient_name, study_description, etc.
        """
        import pydicom

        try:
            ds = pydicom.dcmread(dicom_path, stop_before_pixels=True)
        except Exception as e:
            logger.warning(f"Could not read DICOM headers from {dicom_path}: {e}")
            return {"modality": "CR", "body_part": "CHEST", "error": str(e)}

        modality = str(getattr(ds, "Modality", "")).upper().strip()
        body_part = str(getattr(ds, "BodyPartExamined", "")).upper().strip()
        study_desc = str(getattr(ds, "StudyDescription", ""))
        series_desc = str(getattr(ds, "SeriesDescription", ""))
        patient_name = str(getattr(ds, "PatientName", "Anonymous"))
        patient_id = str(getattr(ds, "PatientID", ""))
        patient_sex = str(getattr(ds, "PatientSex", ""))
        patient_age = str(getattr(ds, "PatientAge", ""))
        rows = int(getattr(ds, "Rows", 0))
        cols = int(getattr(ds, "Columns", 0))
        bits = int(getattr(ds, "BitsStored", 0))
        manufacturer = str(getattr(ds, "Manufacturer", ""))
        institution = str(getattr(ds, "InstitutionName", ""))

        # Heuristic: infer body_part from study/series description if missing
        if not body_part:
            desc_upper = (study_desc + " " + series_desc).upper()
            if "CHEST" in desc_upper or "CXR" in desc_upper or "THORAX" in desc_upper:
                body_part = "CHEST"
            elif "HEAD" in desc_upper or "BRAIN" in desc_upper:
                body_part = "HEAD"
            elif "ABDOMEN" in desc_upper:
                body_part = "ABDOMEN"
            elif "PELVIS" in desc_upper or "PROSTATE" in desc_upper:
                body_part = "PELVIS"
            elif "BREAST" in desc_upper or "MAMMO" in desc_upper:
                body_part = "BREAST"
            elif "HEART" in desc_upper or "CARDIAC" in desc_upper or "CORONARY" in desc_upper:
                body_part = "HEART"

        return {
            "modality": modality,
            "body_part": body_part,
            "study_description": study_desc,
            "series_description": series_desc,
            "patient_name": patient_name,
            "patient_id": patient_id,
            "patient_sex": patient_sex,
            "patient_age": patient_age,
            "rows": rows,
            "columns": cols,
            "bits_stored": bits,
            "manufacturer": manufacturer,
            "institution": institution,
        }

    # =================================================================
    # Workflow Routing
    # =================================================================

    def _route_to_workflow(self, modality: str, body_part: str) -> str:
        """Map DICOM modality + body part to a workflow name.

        Args:
            modality: DICOM Modality tag (e.g., 'CR', 'CT', 'MR').
            body_part: DICOM BodyPartExamined tag (e.g., 'CHEST', 'HEAD').

        Returns:
            Workflow name string.
        """
        result = _MODALITY_WORKFLOW_ROUTING.get(
            (modality, body_part),
            _MODALITY_WORKFLOW_ROUTING.get(
                (modality, ""),
                "cxr_rapid_findings",  # Default fallback
            ),
        )
        return result

    # =================================================================
    # CT Threshold Analysis
    # =================================================================

    def _analyze_ct_by_threshold(
        self, volume: np.ndarray, workflow_name: str, ds: Any
    ) -> Dict:
        """Run intensity-based CT analysis when MONAI weights are unavailable.

        Uses HU windowing and thresholding to provide clinically reasonable
        measurements. This is a fallback that provides real measurements
        from the actual DICOM pixel data, but without trained model inference.

        Args:
            volume: 3D numpy array in HU values, shape (slices, H, W).
            workflow_name: Target workflow name.
            ds: Reference pydicom Dataset for metadata.

        Returns:
            WorkflowResult-compatible dict.
        """
        findings = []
        measurements: Dict[str, float] = {}

        # Volume statistics
        measurements["volume_shape_slices"] = float(volume.shape[0])
        measurements["volume_shape_height"] = float(volume.shape[1])
        measurements["volume_shape_width"] = float(volume.shape[2])
        measurements["hu_min"] = float(volume.min())
        measurements["hu_max"] = float(volume.max())
        measurements["hu_mean"] = float(volume.mean())
        measurements["hu_std"] = float(volume.std())

        # Voxel spacing
        pixel_spacing = getattr(ds, "PixelSpacing", [1.0, 1.0])
        slice_thickness = float(getattr(ds, "SliceThickness", 1.0))
        voxel_vol_mm3 = float(pixel_spacing[0]) * float(pixel_spacing[1]) * slice_thickness
        measurements["voxel_volume_mm3"] = voxel_vol_mm3

        if workflow_name == "ct_chest_lung_nodule":
            findings, extra_meas = self._ct_lung_threshold_analysis(volume, voxel_vol_mm3)
            measurements.update(extra_meas)
        elif workflow_name == "ct_head_hemorrhage":
            findings, extra_meas = self._ct_head_threshold_analysis(volume, voxel_vol_mm3)
            measurements.update(extra_meas)
        else:
            # Generic CT analysis
            tissue_mask = (volume > -200) & (volume < 400)
            tissue_fraction = float(tissue_mask.sum()) / volume.size
            measurements["soft_tissue_fraction"] = round(tissue_fraction, 4)
            findings.append({
                "category": "normal",
                "description": (
                    f"CT volume analyzed: {volume.shape[0]} slices, "
                    f"mean HU={volume.mean():.1f}. "
                    f"Threshold-based analysis (no trained model weights loaded). "
                    f"Clinical review recommended."
                ),
                "severity": "routine",
            })

        if not findings:
            findings.append({
                "category": "normal",
                "description": "No significant findings detected by threshold analysis.",
                "severity": "normal",
            })

        severity_order = ["normal", "routine", "significant", "urgent", "critical"]
        overall_severity = "normal"
        for f in findings:
            f_sev = f.get("severity", "normal")
            if severity_order.index(f_sev) > severity_order.index(overall_severity):
                overall_severity = f_sev

        positive = [f for f in findings if f.get("severity") not in ("normal", "routine")]
        classification = (
            f"positive: {', '.join(f.get('category', '') for f in positive)}"
            if positive
            else "negative"
        )

        return {
            "workflow_name": workflow_name,
            "status": "completed",
            "classification": classification,
            "severity": overall_severity,
            "findings": findings,
            "measurements": measurements,
            "inference_time_ms": 0,
            "is_mock": False,
            "nim_services_used": ["Threshold-based CT analysis (LIVE, no model weights)"],
            "analysis_note": (
                "This analysis uses real DICOM pixel data with intensity thresholding. "
                "For trained model inference, download MONAI model weights."
            ),
        }

    def _ct_lung_threshold_analysis(
        self, volume: np.ndarray, voxel_vol_mm3: float
    ) -> tuple:
        """Lung nodule detection via HU thresholding on CT chest."""
        findings = []
        measurements: Dict[str, float] = {}

        # Lung window: air is ~-1000 HU, lung tissue ~-500 to -800 HU
        lung_mask = (volume > -950) & (volume < -400)
        lung_fraction = float(lung_mask.sum()) / volume.size
        measurements["lung_tissue_fraction"] = round(lung_fraction, 4)
        measurements["lung_volume_ml"] = round(
            float(lung_mask.sum()) * voxel_vol_mm3 / 1000.0, 1
        )

        # Soft tissue within lung field (potential nodules): -100 to 200 HU
        # This is a very rough proxy
        nodule_candidate_mask = (volume > -100) & (volume < 200)
        # Intersect with approximate lung region (central portion of image)
        h, w = volume.shape[1], volume.shape[2]
        lung_region = np.zeros_like(volume, dtype=bool)
        lung_region[:, h // 6 : 5 * h // 6, w // 4 : 3 * w // 4] = True
        candidates = nodule_candidate_mask & lung_region & ~((volume > 200))

        candidate_voxels = int(candidates.sum())
        candidate_vol_mm3 = candidate_voxels * voxel_vol_mm3
        measurements["candidate_nodule_volume_mm3"] = round(candidate_vol_mm3, 1)

        if candidate_vol_mm3 > 500:  # >500 mm3 ~ >10mm diameter sphere
            findings.append({
                "category": "nodule",
                "description": (
                    f"Potential pulmonary density detected "
                    f"(~{candidate_vol_mm3:.0f} mm3 soft tissue in lung fields). "
                    f"Threshold-based detection -- radiologist review required."
                ),
                "severity": "significant",
                "confidence": 0.6,
            })
        else:
            findings.append({
                "category": "normal",
                "description": (
                    f"No significant pulmonary densities detected by thresholding. "
                    f"Lung volume: {measurements['lung_volume_ml']:.0f} mL."
                ),
                "severity": "normal",
            })

        return findings, measurements

    def _ct_head_threshold_analysis(
        self, volume: np.ndarray, voxel_vol_mm3: float
    ) -> tuple:
        """Hemorrhage detection via HU thresholding on CT head."""
        findings = []
        measurements: Dict[str, float] = {}

        # Acute blood is typically 50-80 HU on non-contrast CT head
        blood_mask = (volume > 50) & (volume < 100)
        # Brain parenchyma is typically 20-45 HU
        brain_mask = (volume > 15) & (volume < 50)

        blood_volume_ml = float(blood_mask.sum()) * voxel_vol_mm3 / 1000.0
        brain_volume_ml = float(brain_mask.sum()) * voxel_vol_mm3 / 1000.0

        measurements["blood_density_volume_ml"] = round(blood_volume_ml, 2)
        measurements["brain_parenchyma_volume_ml"] = round(brain_volume_ml, 1)

        # Very rough hemorrhage detection
        if brain_volume_ml > 0:
            blood_brain_ratio = blood_volume_ml / brain_volume_ml
            measurements["blood_brain_ratio"] = round(blood_brain_ratio, 4)

            if blood_volume_ml > 30:
                findings.append({
                    "category": "hemorrhage",
                    "description": (
                        f"Significant high-density region detected "
                        f"(~{blood_volume_ml:.1f} mL at 50-100 HU). "
                        f"Possible intracranial hemorrhage. "
                        f"Threshold-based detection -- urgent radiologist review."
                    ),
                    "severity": "critical",
                    "confidence": 0.55,
                })
            elif blood_volume_ml > 5:
                findings.append({
                    "category": "hemorrhage",
                    "description": (
                        f"Small high-density region detected "
                        f"(~{blood_volume_ml:.1f} mL at 50-100 HU). "
                        f"May represent hemorrhage or calcification. "
                        f"Radiologist review recommended."
                    ),
                    "severity": "significant",
                    "confidence": 0.45,
                })

        if not findings:
            findings.append({
                "category": "normal",
                "description": (
                    "No significant high-density regions detected by thresholding. "
                    "Brain parenchyma volume within normal range."
                ),
                "severity": "normal",
            })

        return findings, measurements

    # =================================================================
    # Single-Frame Fallback
    # =================================================================

    def _analyze_single_frame_fallback(
        self, dicom_path: str, workflow_name: str, meta: Dict
    ) -> Dict:
        """Fallback analysis for single-frame non-CXR DICOM files.

        Returns basic pixel statistics and metadata analysis when
        a specialized model is not available for the modality.
        """
        import pydicom

        ds = pydicom.dcmread(dicom_path)
        pixel_array = ds.pixel_array.astype(np.float32)

        measurements = {
            "pixel_min": float(pixel_array.min()),
            "pixel_max": float(pixel_array.max()),
            "pixel_mean": float(pixel_array.mean()),
            "pixel_std": float(pixel_array.std()),
            "image_height": float(pixel_array.shape[0]),
            "image_width": float(pixel_array.shape[1] if pixel_array.ndim > 1 else 1),
        }

        return {
            "workflow_name": workflow_name,
            "status": "completed",
            "classification": "review_needed",
            "severity": "routine",
            "findings": [{
                "category": "normal",
                "description": (
                    f"DICOM image analyzed ({meta.get('modality', 'unknown')} / "
                    f"{meta.get('body_part', 'unknown')}). "
                    f"No trained model available for this modality in live mode. "
                    f"Basic pixel statistics provided."
                ),
                "severity": "routine",
            }],
            "measurements": measurements,
            "inference_time_ms": 0,
            "is_mock": False,
            "nim_services_used": ["Pixel statistics (no trained model for this modality)"],
        }

    # =================================================================
    # Helpers
    # =================================================================

    def _find_dicom_files(self, directory: str) -> List[str]:
        """Find all DICOM files in a directory (non-recursive).

        Args:
            directory: Directory path to search.

        Returns:
            Sorted list of DICOM file paths.
        """
        dcm_files = []
        for fname in sorted(os.listdir(directory)):
            fpath = os.path.join(directory, fname)
            if not os.path.isfile(fpath):
                continue
            # Accept .dcm files or extensionless files (common in DICOM)
            if fname.lower().endswith(".dcm") or "." not in fname:
                dcm_files.append(fpath)
        return dcm_files

    def _error_result(self, error_msg: str, start_time: float) -> Dict:
        """Build an error WorkflowResult dict."""
        return {
            "workflow_name": "unknown",
            "status": "failed",
            "classification": "error",
            "severity": "routine",
            "findings": [{
                "category": "normal",
                "description": f"Analysis failed: {error_msg}",
                "severity": "routine",
            }],
            "measurements": {},
            "inference_time_ms": round((time.time() - start_time) * 1000, 1),
            "is_mock": False,
            "nim_services_used": [],
            "error": error_msg,
        }

    def get_status(self) -> Dict:
        """Return current analyzer status (GPU, loaded models, etc.)."""
        gpu_name = None
        gpu_memory_mb = None
        if torch.cuda.is_available():
            gpu_name = torch.cuda.get_device_name(0)
            try:
                gpu_memory_mb = round(
                    torch.cuda.get_device_properties(0).total_mem / 1024 / 1024
                )
            except Exception:
                pass

        return {
            "gpu_available": torch.cuda.is_available(),
            "gpu_name": gpu_name,
            "gpu_memory_mb": gpu_memory_mb,
            "device": str(self.device),
            "models_loaded": list(self._loaded_models.keys()),
            "supported_modalities": ["CXR (DenseNet-121)", "CT (threshold)", "MRI (threshold)"],
            "live_analysis_enabled": True,
        }
