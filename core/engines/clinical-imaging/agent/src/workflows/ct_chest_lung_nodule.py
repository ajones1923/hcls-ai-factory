"""CT Chest Lung Nodule Detection & Classification Workflow.

Reference workflow #2: Lung cancer screening nodule analysis.
Uses RetinaNet for detection + SegResNet for volumetric segmentation (MONAI).
Applies ACR Lung-RADS v2022 classification for standardized reporting.

Real inference pipeline (mock_mode=False):
    1. Download MONAI Model Zoo ``lung_nodule_ct_detection`` bundle (RetinaNet)
    2. Load pretrained weights from bundle checkpoint
    3. Run forward pass on a synthetic 3D volume for pipeline validation
    4. Fall back to mock results if download/load fails
"""

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch
from loguru import logger

from src.models import FindingCategory, FindingSeverity, LungRADS, WorkflowResult, WorkflowStatus
from src.workflows.base import BaseImagingWorkflow

# Directory for downloaded MONAI bundles
MONAI_BUNDLE_DIR = os.environ.get("MONAI_BUNDLE_DIR", "/tmp/monai_bundles")
LUNG_NODULE_BUNDLE = "lung_nodule_ct_detection"


# ACR Lung-RADS v2022 thresholds for solid nodules (mm)
LUNG_RADS_THRESHOLDS = {
    "cat_2_max": 6.0,    # Category 2: < 6mm solid
    "cat_3_max": 8.0,    # Category 3: 6-8mm solid
    "cat_4a_max": 15.0,  # Category 4A: 8-15mm solid
    # Category 4B: > 15mm solid
}


def classify_solid_nodule_lung_rads(diameter_mm: float) -> LungRADS:
    """Classify a solid nodule using ACR Lung-RADS v2022 size thresholds.

    Args:
        diameter_mm: Maximum axial diameter of the solid nodule in mm.

    Returns:
        LungRADS category enum value.

    Reference:
        ACR Lung-RADS v2022 Assessment Categories:
        - Category 2 (Benign Appearance): Solid nodule < 6mm
        - Category 3 (Probably Benign): Solid nodule 6-8mm
        - Category 4A (Suspicious): Solid nodule 8-15mm
        - Category 4B (Very Suspicious): Solid nodule >= 15mm

    LIMITATION (documented 2026-11-03): these are the BASELINE-SCREEN thresholds, and size is
    the only input. Lung-RADS v2022 also keys on the exam context, which this function cannot
    see: at follow-up a NEW solid nodule of 4-6mm is category 3 rather than 2, a new 6-8mm
    nodule is 4A rather than 3, and a growing nodule moves up a category. So on a follow-up
    LDCT this function can under-call by one category. The workflow has no prior exam to
    compare against; supplying one is the fix, not changing these thresholds.
    """
    if diameter_mm < LUNG_RADS_THRESHOLDS["cat_2_max"]:
        return LungRADS.CAT_2
    elif diameter_mm < LUNG_RADS_THRESHOLDS["cat_3_max"]:
        return LungRADS.CAT_3
    elif diameter_mm < LUNG_RADS_THRESHOLDS["cat_4a_max"]:
        return LungRADS.CAT_4A
    else:
        return LungRADS.CAT_4B


class LungRADSGradingError(ValueError):
    """Raised when a nodule cannot be graded with the tables this module implements."""


# Nodule densities Lung-RADS v2022 grades, and which of them this module has a table for.
#
# AUDIT FIX (2026-11-03): `classify_solid_nodule_lung_rads` -- a function whose name,
# docstring and published source all say SOLID -- was being handed every nodule, with only the
# choice of WHICH diameter varying by density. Lung-RADS v2022 grades part-solid nodules on
# criteria that differ from the solid table and non-solid (ground-glass) nodules on a different
# scale again, so the solid table was wrong in both directions:
#
#   * a 25 mm part-solid nodule with a 9.3 mm solid component -- the classic invasive
#     adenocarcinoma appearance, and the exact nodule this workflow ships as its demo -- came
#     back 4A, "3-month LDCT or PET/CT", where the published part-solid criteria escalate a
#     solid component of that size immediately;
#   * a 20 mm pure ground-glass nodule, which is indolent and stays on annual screening, came
#     back 4B / CRITICAL / "tissue sampling recommended" -- a biopsy of a lesion that did not
#     need one.
#
# The replacement thresholds are NOT in this file, because the audit could not retrieve them
# from the issuing body: the ACR's own Lung-RADS v2022 PDF redirects to a tenant page and
# serves no document. Guessing them would be the same failure as the defect. So the density is
# now explicit and required, and a density whose published table this module does not hold is
# REFUSED rather than graded with the wrong one. The open question and the document that
# settles it are recorded in CLINICAL_REVIEW.md §6d.
SOLID = "solid"
PART_SOLID = "part_solid"
NON_SOLID = "non_solid"

# Spelling variants seen in this engine's own inputs (the report parser extracts hyphenated
# free-text qualifiers such as "part-solid" and "ground glass"). Normalisation only -- no
# density is inferred from anything other than an explicit statement of it.
_NODULE_TYPE_ALIASES: Dict[str, str] = {
    "solid": SOLID,
    "part_solid": PART_SOLID,
    "partly_solid": PART_SOLID,
    "subsolid": PART_SOLID,
    "sub_solid": PART_SOLID,
    "non_solid": NON_SOLID,
    "nonsolid": NON_SOLID,
    "ground_glass": NON_SOLID,
    "ground_glass_nodule": NON_SOLID,
    "ggn": NON_SOLID,
    "pure_ground_glass": NON_SOLID,
    "ground_glass_opacity": NON_SOLID,
}

# The densities this module can actually grade. Keep this honest: adding a key here without
# adding its published table is how the original defect was written.
_IMPLEMENTED_DENSITY_TABLES = (SOLID,)


def normalize_nodule_type(raw: Any) -> str:
    """Canonicalise a nodule-density label's spelling (case, hyphens, spaces)."""
    text = re.sub(r"[\s\-/]+", "_", str(raw).strip().casefold())
    return re.sub(r"_+", "_", text).strip("_")


def resolve_nodule_type(raw: Any) -> str:
    """Resolve a nodule-density label to one of SOLID / PART_SOLID / NON_SOLID.

    Raises:
        LungRADSGradingError: if the density was not stated, or the label is not one this
            module recognises. Neither is defaulted to "solid": density changes the table.
    """
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        raise LungRADSGradingError(
            "Nodule density was not reported. Lung-RADS v2022 uses a different table for "
            "solid, part-solid and non-solid (ground-glass) nodules, so density cannot be "
            "assumed -- least of all assumed to be solid, which this code used to do."
        )

    key = normalize_nodule_type(raw)
    if key not in _NODULE_TYPE_ALIASES:
        raise LungRADSGradingError(
            f"Unrecognised nodule density {str(raw)!r}. Expected a solid, part-solid or "
            "non-solid (ground-glass) nodule. An unrecognised density is not graded, because "
            "the only way to grade it would be to pick a table at random."
        )
    return _NODULE_TYPE_ALIASES[key]


def classify_nodule_lung_rads(
    nodule_type: Any,
    max_diameter_mm: Optional[float] = None,
    solid_component_diameter_mm: Optional[float] = None,
) -> Tuple[LungRADS, float]:
    """Grade a nodule of any density, or refuse if the published table is not implemented.

    Args:
        nodule_type: The nodule's density, as reported. Required -- see ``resolve_nodule_type``.
        max_diameter_mm: Maximum axial diameter in mm, or None if not measured.
        solid_component_diameter_mm: Solid-component diameter for a part-solid nodule.

    Returns:
        (category, the diameter in mm that the category was read from).

    Raises:
        LungRADSGradingError: if the density is absent, unrecognised, or one whose v2022 table
            this module does not hold; or if the diameter the table needs was not measured.
    """
    density = resolve_nodule_type(nodule_type)

    if density not in _IMPLEMENTED_DENSITY_TABLES:
        raise LungRADSGradingError(
            f"This module has no Lung-RADS v2022 table for a {density.replace('_', '-')} "
            "nodule, so it will not grade one. Applying the solid-nodule table to it mis-calls "
            "in both directions (see the note above this function and CLINICAL_REVIEW.md "
            "§6d). Settle the thresholds against Christensen J et al., 'ACR Lung-RADS v2022: "
            "Assessment Categories and Management Recommendations', Chest 2024;165(3):738-753 "
            "(also J Am Coll Radiol 2024;21(3):473-488), Table 1, then add the table here."
        )

    if max_diameter_mm is None:
        raise LungRADSGradingError(
            "Solid nodule has no measured diameter. Every Lung-RADS size band is a diameter "
            "range, and an absent measurement used to be scored as 0.0 mm -- below every "
            "threshold -- which reported the nodule as category 2, 'benign appearance'."
        )

    return classify_solid_nodule_lung_rads(float(max_diameter_mm)), float(max_diameter_mm)


def lung_rads_to_severity(category: LungRADS) -> FindingSeverity:
    """Map Lung-RADS category to clinical severity."""
    mapping = {
        LungRADS.CAT_1: FindingSeverity.NORMAL,
        LungRADS.CAT_2: FindingSeverity.ROUTINE,
        LungRADS.CAT_3: FindingSeverity.SIGNIFICANT,
        LungRADS.CAT_4A: FindingSeverity.URGENT,
        LungRADS.CAT_4B: FindingSeverity.CRITICAL,
        LungRADS.CAT_4X: FindingSeverity.CRITICAL,
    }
    return mapping.get(category, FindingSeverity.ROUTINE)


def lung_rads_recommendation(category: LungRADS) -> str:
    """Return ACR Lung-RADS v2022 management recommendation."""
    recommendations = {
        LungRADS.CAT_1: "Continue annual LDCT screening.",
        LungRADS.CAT_2: (
            "Continue annual LDCT screening. "
            "Solid nodule < 6mm — benign appearance."
        ),
        LungRADS.CAT_3: (
            "6-month LDCT follow-up recommended. "
            "Solid nodule 6-8mm — probably benign."
        ),
        LungRADS.CAT_4A: (
            "3-month LDCT follow-up or PET/CT recommended. "
            "Solid nodule 8-15mm — suspicious."
        ),
        LungRADS.CAT_4B: (
            "Chest CT with contrast, PET/CT, and/or tissue sampling recommended. "
            "Solid nodule >= 15mm — very suspicious."
        ),
        LungRADS.CAT_4X: (
            "Tissue sampling and/or multidisciplinary consultation recommended. "
            "Additional suspicious features present."
        ),
    }
    return recommendations.get(category, "Clinical correlation recommended.")


class CTChestLungNoduleWorkflow(BaseImagingWorkflow):
    """Lung nodule detection, segmentation, and Lung-RADS classification.

    Pipeline:
        1. Load DICOM CT chest series (low-dose or diagnostic)
        2. Resample to 1mm isotropic, apply lung window (-1000 to 400 HU)
        3. Run RetinaNet nodule detection (bounding boxes + confidence)
        4. Run SegResNet volumetric segmentation on detected nodules
        5. Measure diameter (mm), volume (mm3), density characteristics
        6. Apply ACR Lung-RADS v2022 classification per nodule

    In real mode (mock_mode=False), downloads the MONAI Model Zoo
    ``lung_nodule_ct_detection`` bundle and loads pretrained RetinaNet
    weights trained on the LUNA16 dataset.
    """

    WORKFLOW_NAME: str = "ct_chest_lung_nodule"
    TARGET_LATENCY_SEC: float = 300.0
    MODALITY: str = "ct"
    BODY_REGION: str = "chest"
    MODELS_USED: List[str] = [
        "RetinaNet (MONAI lung_nodule_ct_detection)",
        "SegResNet (MONAI)",
        "VISTA-3D (optional)",
    ]

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        mock_overrides: Optional[Dict] = None,
    ):
        super().__init__(mock_mode=mock_mode, nim_clients=nim_clients, mock_overrides=mock_overrides)
        self._model: Optional[torch.nn.Module] = None
        self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._weights_loaded = False
        self._bundle_name = LUNG_NODULE_BUNDLE

        if not mock_mode:
            self._load_model()

    # ------------------------------------------------------------------
    # Model loading
    # ------------------------------------------------------------------

    def _load_model(self) -> None:
        """Download the MONAI lung_nodule_ct_detection bundle and load the model.

        The bundle contains a RetinaNet architecture trained on LUNA16 for
        3D lung-nodule detection.  If the download or weight-loading fails
        we fall back to a randomly-initialised MONAI RetinaNet so the
        pipeline can still validate end-to-end.
        """
        try:
            self._download_and_load_bundle()
        except Exception as e:
            logger.warning(
                f"Bundle download/load failed ({e}). "
                "Falling back to random-weight RetinaNet."
            )
            self._load_fallback_model()

    def _download_and_load_bundle(self) -> None:
        """Download the bundle via ``monai.bundle.download`` and instantiate the model."""
        import monai.bundle

        bundle_dir = Path(MONAI_BUNDLE_DIR)
        bundle_dir.mkdir(parents=True, exist_ok=True)

        bundle_path = bundle_dir / self._bundle_name

        # Download only if the bundle directory is missing
        if not bundle_path.exists():
            logger.info(
                f"Downloading MONAI bundle '{self._bundle_name}' to {bundle_dir}"
            )
            monai.bundle.download(
                name=self._bundle_name,
                bundle_dir=str(bundle_dir),
            )
            logger.info(f"Bundle downloaded to {bundle_path}")
        else:
            logger.info(f"Bundle already cached at {bundle_path}")

        # --- Try to load via bundle API (ConfigParser) --------------------
        try:
            self._load_via_bundle_config(bundle_path)
            return
        except Exception as cfg_err:
            logger.warning(f"Bundle ConfigParser load failed: {cfg_err}")

        # --- Fallback: manually locate a .pt / .pth checkpoint ------------
        self._load_checkpoint_manually(bundle_path)

    def _load_via_bundle_config(self, bundle_path: Path) -> None:
        """Parse the bundle ``configs/inference.json`` and build the model."""
        from monai.bundle import ConfigParser

        config_path = bundle_path / "configs" / "inference.json"
        if not config_path.exists():
            config_path = bundle_path / "configs" / "inference.yaml"
        if not config_path.exists():
            raise FileNotFoundError(
                f"No inference config found in {bundle_path / 'configs'}"
            )

        parser = ConfigParser()
        parser.read_config(str(config_path))

        # The model key is usually "network_def" or "network"
        for key in ("network_def", "network"):
            try:
                model = parser.get_parsed_content(key)
                break
            except Exception:
                continue
        else:
            raise RuntimeError("Could not resolve model from bundle config")

        # Load checkpoint weights
        ckpt = self._find_checkpoint(bundle_path)
        if ckpt is not None:
            state = torch.load(str(ckpt), map_location=self._device, weights_only=False)
            if isinstance(state, dict) and "state_dict" in state:
                state = state["state_dict"]
            model.load_state_dict(state, strict=False)
            self._weights_loaded = True
            logger.info(f"Loaded pretrained weights from {ckpt.name}")
        else:
            logger.info("No checkpoint found — using config-initialised weights")
            self._weights_loaded = False

        model.to(self._device)
        model.eval()
        self._model = model

        param_count = sum(p.numel() for p in model.parameters())
        logger.info(
            f"Lung-nodule RetinaNet ready: {param_count:,} params, "
            f"device={self._device}, weights_loaded={self._weights_loaded}"
        )

    def _load_checkpoint_manually(self, bundle_path: Path) -> None:
        """Scan the bundle directory for a checkpoint file and load it into a
        fallback MONAI RetinaNet architecture."""
        ckpt = self._find_checkpoint(bundle_path)
        if ckpt is None:
            raise FileNotFoundError(
                f"No .pt/.pth/.ckpt/.ts checkpoint in {bundle_path}"
            )

        logger.info(f"Loading checkpoint directly: {ckpt}")
        state = torch.load(str(ckpt), map_location=self._device, weights_only=False)
        if isinstance(state, dict) and "state_dict" in state:
            state = state["state_dict"]

        # Build a default RetinaNet matching the bundle's architecture
        from monai.apps.detection.networks.retinanet_detector import RetinaNetDetector
        from monai.apps.detection.networks.retinanet_network import resnet_fpn_feature_extractor
        from monai.networks.nets import resnet10

        backbone = resnet10(spatial_dims=3, n_input_channels=1, feed_forward=False)
        feature_extractor = resnet_fpn_feature_extractor(
            backbone=backbone,
            spatial_dims=3,
            pretrained_backbone=False,
            returned_layers=[1, 2],
        )
        # Attempt weight load; may fail on shape mismatch
        try:
            model = RetinaNetDetector(
                network=feature_extractor,
                anchor_generator=None,
                spatial_dims=3,
                num_classes=1,
            )
            model.load_state_dict(state, strict=False)
            self._weights_loaded = True
        except Exception as e:
            logger.warning(f"Weight load mismatch: {e}")
            self._weights_loaded = False

        model.to(self._device)
        model.eval()
        self._model = model

        param_count = sum(p.numel() for p in model.parameters())
        logger.info(
            f"Lung-nodule RetinaNet (manual load) ready: {param_count:,} params, "
            f"weights_loaded={self._weights_loaded}"
        )

    def _load_fallback_model(self) -> None:
        """Create a lightweight MONAI SegResNet with random weights as fallback.

        This gives the pipeline a real PyTorch model to exercise even when the
        full bundle is unavailable.
        """
        from monai.networks.nets import SegResNet

        logger.info("Initialising fallback SegResNet (random weights) for pipeline validation")

        self._model = SegResNet(
            spatial_dims=3,
            in_channels=1,
            out_channels=2,  # background + nodule
            init_filters=8,
            blocks_down=[1, 1, 1],
            blocks_up=[1, 1],
        )
        self._model.to(self._device)
        self._model.eval()
        self._weights_loaded = False

        param_count = sum(p.numel() for p in self._model.parameters())
        logger.info(
            f"Fallback SegResNet ready: {param_count:,} params, device={self._device}"
        )

    @staticmethod
    def _find_checkpoint(bundle_path: Path) -> Optional[Path]:
        """Find the first checkpoint file in a bundle directory tree."""
        models_dir = bundle_path / "models"
        search_dirs = [models_dir, bundle_path] if models_dir.exists() else [bundle_path]
        for d in search_dirs:
            for ext in ("*.pt", "*.pth", "*.ckpt", "*.ts"):
                hits = sorted(d.glob(ext))
                if hits:
                    return hits[0]
        return None

    # ------------------------------------------------------------------
    # Preprocessing
    # ------------------------------------------------------------------

    def preprocess(self, input_path: str) -> Any:
        """Load DICOM CT chest and apply lung-optimized preprocessing.

        Steps:
            - Load DICOM series from input_path
            - Reorient to RAS coordinate system
            - Resample to 1mm isotropic voxel spacing
            - Apply lung window: center=-600 HU, width=1400 HU (range -1000 to 400 HU)
            - Normalize intensity to [0, 1]

        In production, uses monai.transforms:
            LoadImaged, EnsureChannelFirstd, Orientationd, Spacingd,
            ScaleIntensityRanged (a_min=-1000, a_max=400, b_min=0.0, b_max=1.0)
        """
        logger.info(
            f"Preprocessing CT chest from {input_path}: "
            "reorient RAS, resample 1mm iso, lung window -1000 to 400 HU"
        )
        return None

    def infer(self, preprocessed: Any) -> Dict:
        """Run nodule detection (RetinaNet) + segmentation (SegResNet).

        In real mode, creates a synthetic 3D CT volume and runs the MONAI
        model forward pass to validate the pipeline end-to-end.  The
        inference output is combined with mock clinical context (lobe
        assignment, morphology) since full clinical analysis requires
        real DICOM data.

        Falls back to mock results if no model is loaded.
        """
        if self.mock_mode:
            return self._mock_inference()

        if self._model is None:
            logger.warning("No model loaded — falling back to mock inference")
            return self._mock_inference()

        logger.info("Running real MONAI model inference for lung nodule detection")

        try:
            # Create a small synthetic 3D CT volume (lung-window normalised)
            # Shape: [1, 1, 64, 64, 64] — enough to exercise the model
            with torch.no_grad():
                synthetic_volume = torch.randn(
                    1, 1, 64, 64, 64,
                    device=self._device, dtype=torch.float32,
                )

                output = self._model(synthetic_volume)

            # Interpret output depending on model type
            model_info = {
                "model_type": type(self._model).__name__,
                "weights_loaded": self._weights_loaded,
                "bundle": self._bundle_name,
                "device": str(self._device),
                "input_shape": "1x1x64x64x64",
            }

            if isinstance(output, torch.Tensor):
                model_info["output_shape"] = str(list(output.shape))
                model_info["output_dtype"] = str(output.dtype)
            elif isinstance(output, (list, tuple)):
                model_info["output_count"] = len(output)
                for i, o in enumerate(output[:3]):
                    if isinstance(o, torch.Tensor):
                        model_info[f"output_{i}_shape"] = str(list(o.shape))
            elif isinstance(output, dict):
                model_info["output_keys"] = list(output.keys())[:10]

            logger.info(f"Model forward pass succeeded: {model_info}")

            # The forward pass above ran on a SYNTHETIC volume, not on the supplied image, so
            # the clinical findings below are still the canned ones. Say exactly that: a real
            # model ran, and it did not read this patient's pixels.
            mock = self._mock_inference()
            mock["model_info"] = model_info
            mock[self.SYNTHETIC_PASS_KEY] = True
            return mock

        except Exception as e:
            logger.warning(f"Real inference failed ({e}), falling back to mock")
            result = self._mock_inference()
            result["model_inference_error"] = str(e)
            return result

    def _mock_inference(self) -> Dict:
        """Return realistic mock lung nodule detection result.

        Simulates a screening CT finding with two nodules:
            - A 7.2mm solid nodule in the right upper lobe (Lung-RADS 3)
            - A 12.8mm part-solid nodule in the left lower lobe with a 9.3mm solid
              component -- NOT GRADED, because this module holds no Lung-RADS v2022
              part-solid table (see classify_nodule_lung_rads and CLINICAL_REVIEW.md §6d)

        AUDIT FIX (2026-11-03): this docstring called the second nodule "Lung-RADS 4A", which
        is what the code returned by running a part-solid nodule through the solid-nodule
        table. That narrative is what reviewers read instead of running the code, and it is
        what masked the grading defect for a release. The demo payload is deliberately left as
        a part-solid nodule rather than changed to a solid one: the engine's inability to grade
        part-solid disease should be visible in the demo, not hidden by the demo.
        """
        return {
            "nodules_detected": 2,
            "nodules": [
                {
                    "id": "nodule_1",
                    "location": "right upper lobe",
                    "lobe": "RUL",
                    "type": "solid",
                    "max_diameter_mm": 7.2,
                    "min_diameter_mm": 5.8,
                    "mean_diameter_mm": 6.5,
                    "volume_mm3": 195.0,
                    "mean_hu": -12.0,
                    "detection_confidence": 0.94,
                    "centroid_voxel": [128, 200, 310],
                    "bounding_box": [120, 192, 302, 136, 208, 318],
                    "segmentation_mask": None,
                    "morphology": "smooth",
                    "calcification": False,
                    "spiculation": False,
                },
                {
                    "id": "nodule_2",
                    "location": "left lower lobe",
                    "lobe": "LLL",
                    "type": "part-solid",
                    "max_diameter_mm": 12.8,
                    "min_diameter_mm": 10.1,
                    "mean_diameter_mm": 11.5,
                    "volume_mm3": 1102.0,
                    "solid_component_diameter_mm": 9.3,
                    "ground_glass_diameter_mm": 12.8,
                    "mean_hu": 24.0,
                    "detection_confidence": 0.97,
                    "centroid_voxel": [340, 280, 195],
                    "bounding_box": [328, 268, 183, 352, 292, 207],
                    "segmentation_mask": None,
                    "morphology": "lobulated",
                    "calcification": False,
                    "spiculation": True,
                },
            ],
        }

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Classify each nodule via ACR Lung-RADS v2022 and determine overall severity.

        For solid nodules, classification is based on maximum axial diameter.
        For part-solid nodules, the solid component diameter is used.
        Overall severity is driven by the highest-category nodule.
        """
        nodules = inference_result.get("nodules", [])
        nodules_detected = inference_result.get("nodules_detected", len(nodules))

        if nodules_detected == 0 or not nodules:
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.COMPLETED,
                findings=[{
                    "category": FindingCategory.NORMAL.value,
                    "description": "No pulmonary nodules detected",
                    "severity": FindingSeverity.NORMAL.value,
                    "lung_rads": LungRADS.CAT_1.value,
                }],
                measurements={"nodule_count": 0},
                classification=f"Lung-RADS {LungRADS.CAT_1.value}",
                severity=FindingSeverity.NORMAL,
                nim_services_used=self.MODELS_USED,
            )

        findings = []
        overall_severity = FindingSeverity.NORMAL
        highest_lung_rads: Optional[LungRADS] = None
        ungraded: List[str] = []

        for nodule in nodules:
            nodule_id = nodule.get("id", "unknown")
            # AUDIT FIX (2026-11-03): `type` defaulted to "solid" and `max_diameter_mm` to
            # 0.0. Both defaults answered a question the study had not answered: the first
            # chose a grading table, the second put the nodule below every size band.
            nodule_type = nodule.get("type")
            location = nodule.get("location", "unspecified")
            raw_diameter = nodule.get("max_diameter_mm")
            max_diameter = None if raw_diameter is None else float(raw_diameter)

            # Build per-nodule measurement dict
            nodule_measurements: Dict[str, Any] = {
                "max_diameter_mm": max_diameter,
                "volume_mm3": nodule.get("volume_mm3"),
                "detection_confidence": nodule.get("detection_confidence", 0.0),
            }
            if nodule.get("solid_component_diameter_mm") is not None:
                nodule_measurements["solid_component_mm"] = nodule[
                    "solid_component_diameter_mm"
                ]
            if nodule.get("ground_glass_diameter_mm") is not None:
                nodule_measurements["ground_glass_mm"] = nodule["ground_glass_diameter_mm"]

            morphology_notes = []
            if nodule.get("spiculation"):
                morphology_notes.append("spiculated")
            if nodule.get("calcification"):
                morphology_notes.append("calcified")
            if nodule.get("morphology"):
                morphology_notes.append(nodule["morphology"])
            morphology = ", ".join(morphology_notes) if morphology_notes else "unremarkable"

            type_text = str(nodule_type) if nodule_type else "density not reported"

            try:
                lung_rads, _graded_from_mm = classify_nodule_lung_rads(
                    nodule_type,
                    max_diameter_mm=max_diameter,
                    solid_component_diameter_mm=nodule.get("solid_component_diameter_mm"),
                )
            except LungRADSGradingError as e:
                ungraded.append(f"{nodule_id} ({type_text}): {e}")
                findings.append({
                    "category": FindingCategory.INDETERMINATE.value,
                    "description": (
                        f"{type_text} nodule in {location}: Lung-RADS NOT ASSIGNED -- {e}"
                    ),
                    "severity": FindingSeverity.UNKNOWN.value,
                    "nodule_id": nodule_id,
                    "nodule_type": nodule_type,
                    "location": location,
                    "lobe": nodule.get("lobe", ""),
                    "lung_rads": None,
                    "lung_rads_not_assigned_reason": str(e),
                    "recommendation": (
                        "No Lung-RADS category and no screening interval are issued for this "
                        "nodule. Grade it against the published v2022 table for its density "
                        "before acting on this study."
                    ),
                    "measurements": nodule_measurements,
                    "morphology": morphology,
                })
                continue

            nodule_severity = lung_rads_to_severity(lung_rads)
            recommendation = lung_rads_recommendation(lung_rads)

            # Track overall severity (highest wins)
            severity_order = [
                FindingSeverity.NORMAL,
                FindingSeverity.ROUTINE,
                FindingSeverity.SIGNIFICANT,
                FindingSeverity.URGENT,
                FindingSeverity.CRITICAL,
            ]
            if severity_order.index(nodule_severity) > severity_order.index(overall_severity):
                overall_severity = nodule_severity
            if highest_lung_rads is None or (
                severity_order.index(lung_rads_to_severity(lung_rads))
                > severity_order.index(lung_rads_to_severity(highest_lung_rads))
            ):
                highest_lung_rads = lung_rads

            findings.append({
                "category": FindingCategory.NODULE.value,
                "description": (
                    f"{type_text.replace('-', ' ').replace('_', ' ').title()} nodule in "
                    f"{location}, {max_diameter:.1f}mm max diameter, "
                    f"Lung-RADS {lung_rads.value}"
                ),
                "severity": nodule_severity.value,
                "nodule_id": nodule_id,
                "nodule_type": nodule_type,
                "location": location,
                "lobe": nodule.get("lobe", ""),
                "lung_rads": lung_rads.value,
                "recommendation": recommendation,
                "measurements": nodule_measurements,
                "morphology": morphology,
            })

        # Aggregate measurements. Only values that were actually measured are published here:
        # a 0.0 in this dict is read as a measurement of zero, not as an absent one.
        aggregate_measurements: Dict[str, float] = {
            "nodule_count": float(nodules_detected),
        }
        for i, nodule in enumerate(nodules):
            prefix = f"nodule_{i + 1}"
            if nodule.get("max_diameter_mm") is not None:
                aggregate_measurements[f"{prefix}_diameter_mm"] = float(
                    nodule["max_diameter_mm"]
                )
            if nodule.get("volume_mm3") is not None:
                aggregate_measurements[f"{prefix}_volume_mm3"] = float(nodule["volume_mm3"])

        # Add model/bundle metadata
        model_info = inference_result.get("model_info", {})
        if model_info:
            aggregate_measurements["weights_loaded"] = 1.0 if model_info.get("weights_loaded") else 0.0
        if inference_result.get(self.SYNTHETIC_PASS_KEY):
            # A model ran, on synthetic input. Named so it cannot be read as "measured here".
            aggregate_measurements["synthetic_input_forward_pass"] = 1.0

        # A study in which any detected nodule could not be graded is not a completed screen.
        # Reporting the gradeable nodules and letting the rest fall out of the summary is how
        # an ungraded part-solid nodule turned into a reassuring "Lung-RADS 3" study.
        if ungraded:
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.FAILED,
                findings=findings,
                measurements=aggregate_measurements,
                classification=(
                    "Lung-RADS incomplete -- not all detected nodules could be graded: "
                    + " | ".join(ungraded)
                ),
                severity=FindingSeverity.UNKNOWN,
                nim_services_used=self.MODELS_USED,
            )

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=findings,
            measurements=aggregate_measurements,
            classification=f"Lung-RADS {(highest_lung_rads or LungRADS.CAT_1).value}",
            severity=overall_severity,
            nim_services_used=self.MODELS_USED,
        )
