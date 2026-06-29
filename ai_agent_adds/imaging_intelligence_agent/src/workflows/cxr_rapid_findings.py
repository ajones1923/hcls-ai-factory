"""CXR Rapid Findings Multi-Label Classification Workflow.

Reference workflow #3: Emergency chest X-ray triage.
Uses torchxrayvision's CheXpert-trained DenseNet-121 (densenet121-res224-all)
for clinically meaningful multi-label classification of critical findings,
with fallback to MONAI DenseNet-121 + ImageNet transfer if xrv is unavailable.
Optional VILA-M3 vision-language model for natural language report generation.
Target: sub-30-second triage for emergency department workflow.

Real inference pipeline:
    1. Load CXR image (DICOM, PNG, JPEG) via MONAI/PIL
    2. Resize to 224x224, normalize per model requirements
    3. Model forward pass -> sigmoid per-class probabilities
    4. Per-class confidence thresholds -> positive/negative findings
    5. Severity assessment from highest-severity positive finding

Model priority:
    1. torchxrayvision DenseNet (densenet121-res224-all) -- 18 pathology labels
       trained on CheXpert, MIMIC-CXR, NIH, PadChest, and more
    2. MONAI DenseNet-121 + ImageNet backbone (fallback)
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch
from loguru import logger

from src.models import FindingCategory, FindingSeverity, WorkflowResult, WorkflowStatus
from src.workflows.base import BaseImagingWorkflow


# Per-class confidence thresholds (tuned for high specificity in triage)
CXR_CLASS_THRESHOLDS = {
    "pneumothorax": 0.50,
    "consolidation": 0.60,
    "pleural_effusion": 0.55,
    "cardiomegaly": 0.60,
    "fracture": 0.55,
}

# Ordered class labels (must match model output indices for MONAI fallback)
CXR_CLASS_NAMES = [
    "pneumothorax",
    "consolidation",
    "pleural_effusion",
    "cardiomegaly",
    "fracture",
]

# ── torchxrayvision label mapping ──────────────────────────────────────
# The xrv densenet121-res224-all model outputs 18 pathology logits.
# Full label order:
#   0: Atelectasis, 1: Consolidation, 2: Infiltration, 3: Pneumothorax,
#   4: Edema, 5: Emphysema, 6: Fibrosis, 7: Effusion, 8: Pneumonia,
#   9: Pleural_Thickening, 10: Cardiomegaly, 11: Nodule, 12: Mass,
#   13: Hernia, 14: Lung Lesion, 15: Fracture, 16: Lung Opacity,
#   17: Enlarged Cardiomediastinum
#
# We map our 5 target classes to the corresponding xrv output indices:
XRV_CLASS_INDEX_MAP = {
    "pneumothorax": 3,       # xrv label: 'Pneumothorax'
    "consolidation": 1,      # xrv label: 'Consolidation'
    "pleural_effusion": 7,   # xrv label: 'Effusion'
    "cardiomegaly": 10,      # xrv label: 'Cardiomegaly'
    "fracture": 15,          # xrv label: 'Fracture'
}

# Map CXR classes to finding categories
CXR_CLASS_TO_CATEGORY = {
    "pneumothorax": FindingCategory.PNEUMOTHORAX,
    "consolidation": FindingCategory.CONSOLIDATION,
    "pleural_effusion": FindingCategory.EFFUSION,
    "cardiomegaly": FindingCategory.NORMAL,
    "fracture": FindingCategory.FRACTURE,
}

# Per-class severity when positive (base severity before confidence adjustment)
CXR_CLASS_SEVERITY = {
    "pneumothorax": FindingSeverity.CRITICAL,
    "consolidation": FindingSeverity.URGENT,
    "pleural_effusion": FindingSeverity.SIGNIFICANT,
    "cardiomegaly": FindingSeverity.ROUTINE,
    "fracture": FindingSeverity.URGENT,
}

# Per-class clinical descriptions
CXR_CLASS_DESCRIPTIONS = {
    "pneumothorax": (
        "Pneumothorax detected. Evaluate for tension pneumothorax. "
        "Consider chest tube placement if large or symptomatic."
    ),
    "consolidation": (
        "Pulmonary consolidation identified, suggestive of pneumonia or "
        "atelectasis. Clinical correlation with infection markers recommended."
    ),
    "pleural_effusion": (
        "Pleural effusion detected. Assess for laterality, size, and "
        "clinical correlation with cardiac or infectious etiology."
    ),
    "cardiomegaly": (
        "Cardiomegaly identified (cardiothoracic ratio > 0.5). "
        "Echocardiographic correlation recommended."
    ),
    "fracture": (
        "Rib fracture(s) detected. Evaluate for associated pneumothorax "
        "or hemothorax. Pain management and respiratory monitoring advised."
    ),
}

# ImageNet normalization constants (used only for MONAI fallback path)
IMAGENET_MEAN = [0.485, 0.456, 0.406]
IMAGENET_STD = [0.229, 0.224, 0.225]


class CXRRapidFindingsWorkflow(BaseImagingWorkflow):
    """Multi-label chest X-ray classification for emergency triage.

    Pipeline:
        1. Load CXR image (DICOM or PNG)
        2. Resize to 224x224, normalize per model requirements
        3. Run DenseNet-121 multi-label classification
        4. Apply per-class confidence thresholds
        5. Optional: VILA-M3 generates natural language findings description
        6. Determine overall severity from highest-severity positive finding

    Model selection (automatic):
        - **Preferred**: torchxrayvision DenseNet (densenet121-res224-all)
          Trained on real radiology datasets (CheXpert, MIMIC-CXR, NIH, etc.)
          with 18 pathology labels. Produces clinically meaningful predictions.
        - **Fallback**: MONAI DenseNet-121 with ImageNet pretrained backbone.
          Used only when torchxrayvision is not installed.
    """

    WORKFLOW_NAME: str = "cxr_rapid_findings"
    TARGET_LATENCY_SEC: float = 30.0
    MODALITY: str = "cxr"
    BODY_REGION: str = "chest"
    MODELS_USED: List[str] = ["DenseNet-121 (MONAI)", "VILA-M3 (optional)"]

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        checkpoint_path: Optional[str] = None,
    ):
        super().__init__(mock_mode=mock_mode, nim_clients=nim_clients)
        self._model: Optional[torch.nn.Module] = None
        self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._checkpoint_path = checkpoint_path
        self._weights_loaded = False
        self._using_xrv = False  # Track which model backend is active

        if not mock_mode:
            self._load_model()

    def _load_model(self) -> None:
        """Initialize DenseNet-121 model for CXR classification.

        Tries torchxrayvision first (clinically trained weights), then
        falls back to MONAI DenseNet-121 with ImageNet backbone transfer.
        """
        # Priority 1: torchxrayvision (CheXpert-trained, 18 pathology labels)
        if self._checkpoint_path is None:
            try:
                self._load_xrv_model()
                return
            except Exception as e:
                logger.info(
                    f"torchxrayvision not available, falling back to MONAI: {e}"
                )

        # Priority 2: MONAI DenseNet-121 (ImageNet backbone or checkpoint)
        self._load_monai_model()

    def _load_xrv_model(self) -> None:
        """Load torchxrayvision DenseNet trained on real radiology data.

        Model: densenet121-res224-all
        Training data: CheXpert, MIMIC-CXR, NIH ChestX-ray14, PadChest, etc.
        Output: 18 pathology logits (we extract 5 target classes).
        """
        import torchxrayvision as xrv

        logger.info(
            "Loading torchxrayvision DenseNet (densenet121-res224-all) — "
            "CheXpert-trained weights with 18 pathology labels"
        )

        self._model = xrv.models.DenseNet(weights="densenet121-res224-all")
        self._model.to(self._device)
        self._model.eval()
        self._using_xrv = True
        self._weights_loaded = True

        param_count = sum(p.numel() for p in self._model.parameters())
        logger.info(
            f"torchxrayvision DenseNet ready: {param_count:,} params, "
            f"device={self._device}, 18 pathology outputs, "
            f"mapped to {len(XRV_CLASS_INDEX_MAP)} target classes"
        )

    def _load_monai_model(self) -> None:
        """Initialize MONAI DenseNet-121 (fallback when xrv unavailable).

        Loads the MONAI DenseNet-121 architecture (5-class multi-label output).
        If a checkpoint is provided, loads fine-tuned weights. Otherwise
        uses ImageNet backbone transfer or random initialization.
        """
        from monai.networks.nets import DenseNet121

        logger.info(
            "Initializing MONAI DenseNet-121 for CXR 5-class multi-label "
            "classification (fallback path)"
        )

        self._model = DenseNet121(
            spatial_dims=2,
            in_channels=3,
            out_channels=len(CXR_CLASS_NAMES),
        )
        self._using_xrv = False

        # Try loading fine-tuned checkpoint
        if self._checkpoint_path and Path(self._checkpoint_path).exists():
            try:
                state_dict = torch.load(
                    self._checkpoint_path,
                    map_location=self._device,
                    weights_only=True,
                )
                self._model.load_state_dict(state_dict, strict=False)
                self._weights_loaded = True
                logger.info(f"Loaded CXR fine-tuned weights from {self._checkpoint_path}")
            except Exception as e:
                logger.warning(f"Could not load checkpoint {self._checkpoint_path}: {e}")
                self._weights_loaded = False
        else:
            # Try torchvision ImageNet pretrained backbone
            try:
                self._load_imagenet_backbone()
            except Exception as e:
                logger.info(f"No pretrained backbone available: {e}")
                self._weights_loaded = False

        self._model.to(self._device)
        self._model.eval()

        param_count = sum(p.numel() for p in self._model.parameters())
        logger.info(
            f"MONAI DenseNet-121 ready: {param_count:,} params, device={self._device}, "
            f"weights_loaded={self._weights_loaded}"
        )

    def _load_imagenet_backbone(self) -> None:
        """Load ImageNet-pretrained DenseNet-121 features from torchvision.

        Transfers the pretrained feature extractor weights but keeps the
        classifier head randomly initialized (5 CXR classes != 1000 ImageNet).
        This is standard practice for medical imaging transfer learning.
        """
        try:
            from torchvision.models import densenet121, DenseNet121_Weights

            pretrained = densenet121(weights=DenseNet121_Weights.IMAGENET1K_V1)

            # Copy feature extractor weights (everything except final classifier)
            pretrained_dict = pretrained.state_dict()
            model_dict = self._model.state_dict()

            # Filter: only copy layers that match in name AND shape
            transfer_dict = {
                k: v for k, v in pretrained_dict.items()
                if k in model_dict and v.shape == model_dict[k].shape
            }

            model_dict.update(transfer_dict)
            self._model.load_state_dict(model_dict)

            self._weights_loaded = True
            logger.info(
                f"Loaded ImageNet pretrained backbone: "
                f"{len(transfer_dict)}/{len(model_dict)} layers transferred"
            )
        except ImportError:
            logger.info("torchvision not available — using random initialization")
            self._weights_loaded = False

    def preprocess(self, input_path: str) -> torch.Tensor:
        """Load chest X-ray and preprocess for model inference.

        Dispatches to the appropriate preprocessing pipeline based on
        which model backend is active:
            - xrv: grayscale [1, 1, 224, 224] with xrv.datasets.normalize
            - MONAI: RGB [1, 3, 224, 224] with ImageNet normalization
        """
        input_path = str(input_path)
        logger.info(f"Preprocessing CXR: {input_path}")

        ext = Path(input_path).suffix.lower()

        if ext == ".dcm":
            image = self._load_dicom(input_path)
        elif ext in (".nii", ".gz"):
            image = self._load_nifti(input_path)
        else:
            image = self._load_standard_image(input_path)

        if self._using_xrv:
            return self._preprocess_xrv(image)
        else:
            return self._preprocess_monai(image)

    def _preprocess_xrv(self, image: np.ndarray) -> torch.Tensor:
        """Preprocess image for torchxrayvision DenseNet.

        torchxrayvision expects:
            - Grayscale image normalized to [0, 255] range
            - xrv.datasets.normalize(image, 255) applied
            - Input shape: (batch, 1, 224, 224)
            - The model handles channel replication internally
        """
        import torchxrayvision as xrv

        # Ensure 2D grayscale
        if image.ndim == 3:
            # Take first channel or average
            if image.shape[-1] == 1:
                image = image[..., 0]
            else:
                image = np.mean(image, axis=-1)

        image = image.astype(np.float32)

        # Scale to [0, 255] range for xrv normalization
        if image.max() <= 1.0 and image.max() > 0:
            image = image * 255.0
        elif image.max() > 255.0:
            image = (image / image.max()) * 255.0

        # Resize to 224x224
        image = self._resize(image, target_size=(224, 224))

        # Apply torchxrayvision normalization (centers to [-1024, 1024] range)
        image = xrv.datasets.normalize(image, maxval=255)

        # Convert to tensor: (1, 1, 224, 224)
        tensor = torch.from_numpy(image).float().unsqueeze(0).unsqueeze(0)
        tensor = tensor.to(self._device)

        logger.info(
            f"Preprocessed CXR for xrv: {tensor.shape}, "
            f"range=[{tensor.min():.2f}, {tensor.max():.2f}], "
            f"device={tensor.device}"
        )
        return tensor

    def _preprocess_monai(self, image: np.ndarray) -> torch.Tensor:
        """Preprocess image for MONAI DenseNet-121 (fallback path).

        Pipeline:
            1. Convert to float32 [0, 1] range
            2. Resize to 224x224
            3. Replicate to 3 channels if grayscale
            4. Normalize with ImageNet statistics
            5. Return as torch.Tensor [1, 3, 224, 224]
        """
        # Ensure float32 in [0, 1]
        image = image.astype(np.float32)
        if image.max() > 1.0:
            image = image / image.max() if image.max() > 0 else image

        # Resize to 224x224
        image = self._resize(image, target_size=(224, 224))

        # Ensure 3 channels (H, W, C)
        if image.ndim == 2:
            image = np.stack([image] * 3, axis=-1)
        elif image.ndim == 3 and image.shape[-1] == 1:
            image = np.repeat(image, 3, axis=-1)
        elif image.ndim == 3 and image.shape[-1] > 3:
            image = image[..., :3]

        # ImageNet normalization per channel
        for c in range(3):
            image[..., c] = (image[..., c] - IMAGENET_MEAN[c]) / IMAGENET_STD[c]

        # Convert to tensor: [1, 3, 224, 224]
        tensor = torch.from_numpy(image).permute(2, 0, 1).unsqueeze(0)
        tensor = tensor.to(self._device)

        logger.info(f"Preprocessed CXR for MONAI: {tensor.shape}, device={tensor.device}")
        return tensor

    def _load_dicom(self, path: str) -> np.ndarray:
        """Load a DICOM file to numpy array."""
        try:
            import pydicom
            ds = pydicom.dcmread(path)
            pixel_array = ds.pixel_array.astype(np.float32)
            # Apply windowing if available
            if hasattr(ds, "WindowCenter") and hasattr(ds, "WindowWidth"):
                center = float(ds.WindowCenter) if not isinstance(ds.WindowCenter, float) else ds.WindowCenter
                width = float(ds.WindowWidth) if not isinstance(ds.WindowWidth, float) else ds.WindowWidth
                lower = center - width / 2
                upper = center + width / 2
                pixel_array = np.clip(pixel_array, lower, upper)
                pixel_array = (pixel_array - lower) / (upper - lower)
            return pixel_array
        except Exception as e:
            logger.warning(f"DICOM load failed, falling back to PIL: {e}")
            return self._load_standard_image(path)

    def _load_nifti(self, path: str) -> np.ndarray:
        """Load a NIfTI file (take middle slice if 3D)."""
        try:
            import nibabel as nib
            img = nib.load(path)
            data = img.get_fdata().astype(np.float32)
            if data.ndim == 3:
                mid = data.shape[2] // 2
                data = data[:, :, mid]
            return data
        except Exception as e:
            logger.warning(f"NIfTI load failed: {e}")
            raise

    def _load_standard_image(self, path: str) -> np.ndarray:
        """Load PNG/JPEG/other standard format via PIL."""
        from PIL import Image
        img = Image.open(path).convert("L")  # Convert to grayscale
        return np.array(img, dtype=np.float32)

    def _resize(self, image: np.ndarray, target_size: tuple = (224, 224)) -> np.ndarray:
        """Resize image to target size using PIL (high quality)."""
        from PIL import Image

        if image.ndim == 3:
            # Multi-channel: resize each channel
            resized_channels = []
            for c in range(image.shape[-1]):
                pil_img = Image.fromarray(image[..., c].astype(np.float32), mode="F")
                pil_img = pil_img.resize(target_size, Image.BILINEAR)
                resized_channels.append(np.array(pil_img))
            return np.stack(resized_channels, axis=-1)
        else:
            pil_img = Image.fromarray(image.astype(np.float32), mode="F")
            pil_img = pil_img.resize(target_size, Image.BILINEAR)
            return np.array(pil_img)

    def infer(self, preprocessed: torch.Tensor) -> Dict:
        """Run model forward pass on preprocessed CXR tensor.

        Dispatches to the appropriate inference path:
            - xrv: extract 5 target classes from 18-class output
            - MONAI: direct 5-class output
        Returns per-class sigmoid probabilities for the 5 target findings.
        """
        if self.mock_mode:
            return self._mock_inference()

        if self._model is None:
            raise RuntimeError("Model not loaded. Initialize with mock_mode=False.")

        if self._using_xrv:
            return self._infer_xrv(preprocessed)
        else:
            return self._infer_monai(preprocessed)

    def _infer_xrv(self, preprocessed: torch.Tensor) -> Dict:
        """Run torchxrayvision DenseNet forward pass.

        The model outputs 18 pathology logits. We extract the 5 target
        classes using XRV_CLASS_INDEX_MAP and apply sigmoid activation.
        """
        logger.info(
            "Running torchxrayvision DenseNet classification "
            "(densenet121-res224-all, 18 pathology labels)"
        )

        with torch.no_grad():
            logits = self._model(preprocessed)  # shape: (1, 18)
            all_probs = torch.sigmoid(logits).cpu().numpy()[0]

        # Extract only the 5 target classes from the 18-class output
        class_probs = {}
        for class_name in CXR_CLASS_NAMES:
            xrv_idx = XRV_CLASS_INDEX_MAP[class_name]
            class_probs[class_name] = round(float(all_probs[xrv_idx]), 4)

        logger.info(f"CXR class probabilities (xrv): {class_probs}")

        # Log all 18 pathology scores at debug level for full transparency
        xrv_all_labels = [
            "Atelectasis", "Consolidation", "Infiltration", "Pneumothorax",
            "Edema", "Emphysema", "Fibrosis", "Effusion", "Pneumonia",
            "Pleural_Thickening", "Cardiomegaly", "Nodule", "Mass",
            "Hernia", "Lung Lesion", "Fracture", "Lung Opacity",
            "Enlarged Cardiomediastinum",
        ]
        full_scores = {
            label: round(float(prob), 4)
            for label, prob in zip(xrv_all_labels, all_probs)
        }
        logger.debug(f"Full xrv 18-class scores: {full_scores}")

        return {
            "class_probabilities": class_probs,
            "model_name": "DenseNet-121 (torchxrayvision, densenet121-res224-all)",
            "input_size": "224x224",
            "weights_loaded": self._weights_loaded,
            "xrv_all_scores": full_scores,
        }

    def _infer_monai(self, preprocessed: torch.Tensor) -> Dict:
        """Run MONAI DenseNet-121 forward pass (fallback path)."""
        logger.info("Running MONAI DenseNet-121 CXR classification inference (fallback)")

        with torch.no_grad():
            logits = self._model(preprocessed)
            probabilities = torch.sigmoid(logits).cpu().numpy()[0]

        class_probs = {
            name: round(float(prob), 4)
            for name, prob in zip(CXR_CLASS_NAMES, probabilities)
        }

        logger.info(f"CXR class probabilities (MONAI): {class_probs}")

        return {
            "class_probabilities": class_probs,
            "model_name": "DenseNet-121 (MONAI, ImageNet backbone)",
            "input_size": "224x224",
            "weights_loaded": self._weights_loaded,
        }

    def _mock_inference(self) -> Dict:
        """Return realistic mock multi-label CXR classification result.

        Simulates a case with consolidation (likely pneumonia) and small
        pleural effusion -- a common emergency department presentation.
        """
        return {
            "class_probabilities": {
                "pneumothorax": 0.08,
                "consolidation": 0.87,
                "pleural_effusion": 0.72,
                "cardiomegaly": 0.31,
                "fracture": 0.12,
            },
            "model_name": "DenseNet-121 (torchxrayvision, densenet121-res224-all)",
            "input_size": "224x224",
            "weights_loaded": False,
        }

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Apply per-class thresholds and determine overall severity.

        For each class, the sigmoid probability is compared against the
        class-specific threshold. Positive findings are added to the
        findings list. Overall severity is the maximum severity among
        all positive findings.
        """
        class_probs = inference_result.get("class_probabilities", {})
        weights_loaded = inference_result.get("weights_loaded", False)

        if not class_probs:
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.COMPLETED,
                findings=[{
                    "category": FindingCategory.NORMAL.value,
                    "description": "No classification output available",
                    "severity": FindingSeverity.ROUTINE.value,
                }],
                measurements={},
                classification="incomplete",
                severity=FindingSeverity.ROUTINE,
                nim_services_used=self.MODELS_USED,
            )

        findings = []
        positive_classes = []
        overall_severity = FindingSeverity.NORMAL

        severity_order = [
            FindingSeverity.NORMAL,
            FindingSeverity.ROUTINE,
            FindingSeverity.SIGNIFICANT,
            FindingSeverity.URGENT,
            FindingSeverity.CRITICAL,
        ]

        # Evaluate each class against its threshold
        for class_name, probability in class_probs.items():
            threshold = CXR_CLASS_THRESHOLDS.get(class_name, 0.50)
            is_positive = probability >= threshold

            if is_positive:
                positive_classes.append(class_name)
                class_severity = CXR_CLASS_SEVERITY.get(class_name, FindingSeverity.ROUTINE)
                class_category = CXR_CLASS_TO_CATEGORY.get(class_name, FindingCategory.NORMAL)
                class_description = CXR_CLASS_DESCRIPTIONS.get(class_name, "Finding detected.")

                findings.append({
                    "category": class_category.value,
                    "description": class_description,
                    "severity": class_severity.value,
                    "class_name": class_name,
                    "confidence": round(probability, 3),
                    "threshold": threshold,
                    "above_threshold": True,
                })

                if severity_order.index(class_severity) > severity_order.index(overall_severity):
                    overall_severity = class_severity

        # If no positive findings, report as normal
        if not findings:
            findings.append({
                "category": FindingCategory.NORMAL.value,
                "description": (
                    "No significant acute findings. All classification scores "
                    "below clinical thresholds."
                ),
                "severity": FindingSeverity.NORMAL.value,
            })

        # Build measurements from all class probabilities
        measurements = {}
        for class_name, probability in class_probs.items():
            measurements[f"{class_name}_confidence"] = round(probability, 4)

        measurements["weights_loaded"] = 1.0 if weights_loaded else 0.0
        measurements["using_xrv"] = 1.0 if self._using_xrv else 0.0

        # Classification summary
        if positive_classes:
            classification = f"positive: {', '.join(sorted(positive_classes))}"
        else:
            classification = "negative"

        # Determine which models were used
        model_name = inference_result.get("model_name", "DenseNet-121")
        models_used = [model_name, "VILA-M3 (optional)"]

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=findings,
            measurements=measurements,
            classification=classification,
            severity=overall_severity,
            nim_services_used=models_used,
        )
