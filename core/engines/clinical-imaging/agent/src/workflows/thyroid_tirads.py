"""Thyroid TI-RADS Classification Workflow.

Clinical workflow for thyroid ultrasound with ACR TI-RADS scoring.
Cross-modal trigger: TI-RADS 4+ -> BRAF V600E, RAS, RET/PTC genomic queries.

Scoring system: ACR TI-RADS (2017)
  Points system based on 5 features:
  - Composition: cystic(0), spongiform(0), mixed(1), solid(2)
  - Echogenicity: anechoic(0), hyper/isoechoic(1), hypoechoic(2), very hypoechoic(3)
  - Shape: wider-than-tall(0), taller-than-wide(3)
  - Margin: smooth(0), ill-defined(0), lobulated/irregular(2), extra-thyroidal extension(3)
  - Echogenic foci: none(0), comet-tail(0), macrocalcifications(1), peripheral(2), punctate(3)

  Total points -> TR level:
  TR1 (0 pts): Benign -- no FNA
  TR2 (2 pts): Not suspicious -- no FNA
  TR3 (3 pts): Mildly suspicious -- FNA if >= 2.5cm, follow if >= 1.5cm
  TR4 (4-6 pts): Moderately suspicious -- FNA if >= 1.5cm, follow if >= 1.0cm
  TR5 (7+ pts): Highly suspicious -- FNA if >= 1.0cm, follow if >= 0.5cm

Target latency: < 3 minutes
Models: MONAI thyroid segmentation

Author: Adam Jones
Date: April 2026
"""

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import (
    FindingCategory,
    FindingSeverity,
    TIRADS,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows.base import BaseImagingWorkflow


# =====================================================================
# TI-RADS point values per feature category
# =====================================================================

COMPOSITION_POINTS = {
    "cystic": 0,
    "spongiform": 0,
    "mixed_cystic_solid": 1,
    "solid": 2,
}

ECHOGENICITY_POINTS = {
    "anechoic": 0,
    "hyperechoic": 1,
    "isoechoic": 1,
    "hypoechoic": 2,
    "very_hypoechoic": 3,
}

SHAPE_POINTS = {
    "wider_than_tall": 0,
    "taller_than_wide": 3,
}

MARGIN_POINTS = {
    "smooth": 0,
    "ill_defined": 0,
    "lobulated": 2,
    "irregular": 2,
    "extra_thyroidal_extension": 3,
}

ECHOGENIC_FOCI_POINTS = {
    "none": 0,
    "large_comet_tail": 0,
    "macrocalcifications": 1,
    "peripheral_calcifications": 2,
    "punctate_echogenic_foci": 3,
}


# =====================================================================
# TI-RADS helper functions
# =====================================================================


# The five ACR feature tables, keyed by the nodule field each one scores. Collected so the
# scorer iterates them instead of repeating a lookup per feature -- the repetition is what let
# five identical `.get(..., 0)` defaults go unnoticed.
FEATURE_POINT_TABLES: Dict[str, Dict[str, int]] = {
    "composition": COMPOSITION_POINTS,
    "echogenicity": ECHOGENICITY_POINTS,
    "shape": SHAPE_POINTS,
    "margin": MARGIN_POINTS,
    "echogenic_foci": ECHOGENIC_FOCI_POINTS,
}


class TIRADSFeatureError(ValueError):
    """Raised when a nodule's features cannot be scored against the ACR tables."""


def normalize_feature_label(value: Any) -> str:
    """Canonicalise a feature label's spelling so a lookup miss means a real mismatch.

    Casefolds, trims, and maps hyphens, spaces and slashes to the underscore the tables use,
    so ``"Taller-Than-Wide"``, ``"taller than wide"`` and ``"taller_than_wide"`` are one key.
    This normalises SPELLING only -- it invents no synonyms. A label that still does not match
    after this is an unrecognised finding and is reported as one.
    """
    text = str(value).strip().casefold()
    text = re.sub(r"[\s\-/]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


@dataclass
class TIRADSFeatureScore:
    """The result of scoring one nodule's five ACR features.

    ``points`` is the total of the features that COULD be scored, so it is only a TI-RADS
    total when both ``missing`` and ``unrecognised`` are empty. Check ``is_complete`` before
    mapping it to a TR level.
    """
    points: int = 0
    scored: Dict[str, int] = field(default_factory=dict)
    missing: List[str] = field(default_factory=list)
    unrecognised: Dict[str, str] = field(default_factory=dict)

    @property
    def is_complete(self) -> bool:
        return not self.missing and not self.unrecognised

    def describe_problems(self) -> str:
        """One human-readable clause naming every feature that could not be scored."""
        parts: List[str] = []
        if self.missing:
            parts.append("not reported: " + ", ".join(sorted(self.missing)))
        if self.unrecognised:
            parts.append(
                "not in the ACR tables: "
                + ", ".join(f"{k}={v!r}" for k, v in sorted(self.unrecognised.items()))
            )
        return "; ".join(parts)


def score_tirads_features(nodule: Dict) -> TIRADSFeatureScore:
    """Score a nodule's five ACR TI-RADS features, reporting what could not be scored.

    AUDIT FIX (2026-11-03): this logic used to be five ``dict.get(feature, 0)`` calls. Zero is
    the BENIGN value in all five ACR tables, so an absent key, a vendor synonym, a different
    capitalisation, or a space where the table has an underscore all scored as the single most
    reassuring finding available, silently. A 13-point TR5 nodule written as ``"Solid"``,
    ``"very hypoechoic"``, ``"taller-than-wide"``, ``"Irregular"``, ``"punctate echogenic
    foci"`` came back as 0 points, TR1, "Benign. No FNA indicated. No follow-up needed." --
    and the report still carried a TI-RADS level, so it looked like the system had worked.
    The engine's own report parser already extracts hyphenated free-text qualifiers, so
    non-canonical spellings are this engine's normal input, not a hypothetical.

    Spelling is normalised first (see ``normalize_feature_label``); anything still unmatched is
    returned in ``unrecognised`` and contributes NO points, rather than contributing zero as
    though the finding were benign.

    Args:
        nodule: Dict with keys composition, echogenicity, shape, margin, echogenic_foci.

    Returns:
        TIRADSFeatureScore. Only trust ``.points`` when ``.is_complete`` is True.
    """
    score = TIRADSFeatureScore()

    for feature, table in FEATURE_POINT_TABLES.items():
        raw = nodule.get(feature)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            score.missing.append(feature)
            continue

        key = normalize_feature_label(raw)
        if key not in table:
            score.unrecognised[feature] = str(raw)
            continue

        value = table[key]
        score.scored[feature] = value
        score.points += value

    return score


def calculate_tirads_points(nodule: Dict) -> int:
    """Calculate total ACR TI-RADS points from nodule features.

    Args:
        nodule: Dict with keys composition, echogenicity, shape,
                margin, echogenic_foci.

    Returns:
        Total TI-RADS points (0-14 theoretical max).

    Raises:
        TIRADSFeatureError: if any of the five features is absent or carries a label that is
            not in its ACR table. See ``score_tirads_features`` for why this refuses rather
            than scoring the unknown feature as benign. Callers that need to carry on and
            report an indeterminate result should use ``score_tirads_features`` directly.
    """
    score = score_tirads_features(nodule)
    if not score.is_complete:
        raise TIRADSFeatureError(
            f"Cannot score TI-RADS points -- {score.describe_problems()}. "
            "An unscored feature is not a benign feature: zero points is the most reassuring "
            "value in every ACR table, so it cannot be used as a default."
        )
    return score.points


def points_to_tirads(points: int) -> TIRADS:
    """Map total TI-RADS points to TR level.

    ACR TI-RADS (2017):
        0 pts:   TR1 (Benign)
        2 pts:   TR2 (Not suspicious)
        3 pts:   TR3 (Mildly suspicious)
        4-6 pts: TR4 (Moderately suspicious)
        7+ pts:  TR5 (Highly suspicious)

    Args:
        points: Total TI-RADS points.

    Returns:
        TIRADS enum value.
    """
    if points <= 0:
        return TIRADS.TR1
    elif points <= 2:
        # NOTE (2026-11-03): the published chart names 0, 2, 3, 4-6 and 7+ and does not say
        # what a total of exactly 1 is. One point is reachable (mixed cystic/solid composition
        # with no other suspicious feature), and it is treated as TR2 here, which is the
        # conventional reading -- TR2 is "not suspicious" and neither 1 nor 2 points warrants
        # FNA or follow-up at any size, so the management is identical either way. Recorded in
        # CLINICAL_REVIEW.md rather than left as a silent choice.
        return TIRADS.TR2
    elif points == 3:
        return TIRADS.TR3
    elif points <= 6:
        return TIRADS.TR4
    else:
        return TIRADS.TR5


# TR levels in published order, for picking the highest nodule in a gland.
#
# AUDIT FIX (2026-11-03): the workflow used to track the highest TR level by SEVERITY index,
# and `tirads_to_severity` maps TR1 and TR2 both to NORMAL. A gland whose only nodule was
# 2 points therefore reported "TI-RADS TR1" in the headline while the finding underneath said
# TR2. Management is identical at TR1 and TR2, so this was a reporting-integrity defect rather
# than a management error -- but TR1 means "benign" and TR2 means "not suspicious", which are
# different statements about a real nodule, and a headline contradicting its own body is not
# a report anyone should sign.
TIRADS_ORDER: Dict[TIRADS, int] = {
    TIRADS.TR1: 1,
    TIRADS.TR2: 2,
    TIRADS.TR3: 3,
    TIRADS.TR4: 4,
    TIRADS.TR5: 5,
}


def tirads_to_severity(tr_level: TIRADS) -> FindingSeverity:
    """Map TI-RADS level to clinical severity.

    Args:
        tr_level: ACR TI-RADS level (TR1-TR5).

    Returns:
        FindingSeverity corresponding to the TR level.
    """
    mapping = {
        TIRADS.TR1: FindingSeverity.NORMAL,
        TIRADS.TR2: FindingSeverity.NORMAL,
        TIRADS.TR3: FindingSeverity.SIGNIFICANT,
        TIRADS.TR4: FindingSeverity.URGENT,
        TIRADS.TR5: FindingSeverity.CRITICAL,
    }
    return mapping.get(tr_level, FindingSeverity.ROUTINE)


_TR_LABELS: Dict[TIRADS, str] = {
    TIRADS.TR3: "TR3 (Mildly suspicious)",
    TIRADS.TR4: "TR4 (Moderately suspicious)",
    TIRADS.TR5: "TR5 (Highly suspicious)",
}


def tirads_fna_recommendation(
    tr_level: TIRADS,
    max_diameter_mm: Optional[float],
) -> str:
    """Determine FNA recommendation based on TR level and nodule size.

    ACR TI-RADS size thresholds for FNA:
        TR1: No FNA at any size
        TR2: No FNA at any size
        TR3: FNA if >= 25mm, follow-up if >= 15mm
        TR4: FNA if >= 15mm, follow-up if >= 10mm
        TR5: FNA if >= 10mm, follow-up if >= 5mm

    AUDIT FIX (2026-11-03): ``max_diameter_mm`` was non-optional and the workflow supplied
    ``0.0`` when no measurement had come through. Zero sits below every threshold above, so a
    nodule correctly scored TR5 at 13 points was reported as "TR5 (Highly suspicious). No FNA
    needed at current size (< 0.5cm). Optional follow-up ultrasound." The report asserted a
    size it never had, and that assertion is the one that cancels the biopsy. Size is now
    Optional and an absent size produces an explicit refusal naming the missing measurement.
    TR1 and TR2 still answer, because at those levels the published management is the same at
    every size and nothing turns on the number.

    Args:
        tr_level: ACR TI-RADS level.
        max_diameter_mm: Maximum nodule diameter in mm, or None if it was not measured.

    Returns:
        Management recommendation string.
    """
    if tr_level == TIRADS.TR1:
        return "TR1 (Benign). No FNA indicated. No follow-up needed."

    if tr_level == TIRADS.TR2:
        return "TR2 (Not suspicious). No FNA indicated. No follow-up needed."

    if max_diameter_mm is None:
        label = _TR_LABELS.get(tr_level, str(tr_level.value))
        return (
            f"{label}. Maximum diameter not measured, and every FNA and follow-up threshold "
            "at this level is size-dependent, so management CANNOT be determined. Measure the "
            "nodule and re-run. This is not a recommendation against FNA."
        )

    if tr_level == TIRADS.TR3:
        if max_diameter_mm >= 25.0:
            return (
                "TR3 (Mildly suspicious). FNA recommended (nodule >= 2.5cm). "
                "Consider ultrasound-guided fine needle aspiration."
            )
        elif max_diameter_mm >= 15.0:
            return (
                "TR3 (Mildly suspicious). Follow-up ultrasound recommended "
                "at 1, 3, and 5 years (nodule >= 1.5cm but < 2.5cm)."
            )
        else:
            return (
                "TR3 (Mildly suspicious). No FNA or follow-up needed "
                "at current size (< 1.5cm)."
            )

    if tr_level == TIRADS.TR4:
        if max_diameter_mm >= 15.0:
            return (
                "TR4 (Moderately suspicious). FNA recommended (nodule >= 1.5cm). "
                "Ultrasound-guided fine needle aspiration indicated."
            )
        elif max_diameter_mm >= 10.0:
            return (
                "TR4 (Moderately suspicious). Follow-up ultrasound recommended "
                "at 1, 2, 3, and 5 years (nodule >= 1.0cm but < 1.5cm)."
            )
        else:
            return (
                "TR4 (Moderately suspicious). No FNA or follow-up needed "
                "at current size (< 1.0cm)."
            )

    if tr_level == TIRADS.TR5:
        if max_diameter_mm >= 10.0:
            return (
                "TR5 (Highly suspicious). FNA recommended (nodule >= 1.0cm). "
                "Urgent ultrasound-guided fine needle aspiration indicated."
            )
        elif max_diameter_mm >= 5.0:
            return (
                "TR5 (Highly suspicious). Follow-up ultrasound recommended "
                "annually for up to 5 years (nodule >= 0.5cm but < 1.0cm)."
            )
        else:
            return (
                "TR5 (Highly suspicious). No FNA needed at current size "
                "(< 0.5cm). Optional follow-up ultrasound."
            )

    return "Clinical correlation recommended."


# =====================================================================
# ThyroidTIRADSWorkflow
# =====================================================================


class ThyroidTIRADSWorkflow(BaseImagingWorkflow):
    """Thyroid nodule TI-RADS classification workflow.

    Pipeline:
        1. Load thyroid ultrasound DICOM series
        2. Segment thyroid gland and identify nodules (MONAI)
        3. Extract nodule features (composition, echogenicity, shape, margin, foci)
        4. Calculate TI-RADS points per nodule
        5. Map points to TR level (TR1-TR5)
        6. Determine FNA recommendation based on TR level + nodule size
        7. Trigger cross-modal genomic query if TR >= 4

    Cross-modal trigger:
        TI-RADS 4+ fires BRAF V600E, RAS, RET/PTC genomic queries
        via CrossModalTrigger to enrich findings with thyroid cancer
        molecular marker context.
    """

    WORKFLOW_NAME: str = "thyroid_tirads"
    TARGET_LATENCY_SEC: float = 180.0
    MODALITY: str = "ultrasound"
    BODY_REGION: str = "neck"
    MODELS_USED: List[str] = ["monai_thyroid_segmentation"]

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        mock_overrides: Optional[Dict] = None,
    ):
        super().__init__(
            mock_mode=mock_mode,
            nim_clients=nim_clients,
            mock_overrides=mock_overrides,
        )

    # ------------------------------------------------------------------
    # Preprocessing
    # ------------------------------------------------------------------

    def preprocess(self, input_path: str) -> Any:
        """Load thyroid ultrasound DICOM and apply preprocessing.

        Steps:
            - Load DICOM ultrasound series from input_path
            - Apply speckle reduction filtering
            - Normalize intensity to [0, 1]
            - Segment thyroid gland boundary

        In production, uses monai.transforms:
            LoadImaged, EnsureChannelFirstd, ScaleIntensityRanged,
            ultrasound-specific speckle noise reduction.
        """
        logger.info(
            f"Preprocessing thyroid ultrasound from {input_path}: "
            "speckle reduction, intensity normalization, gland segmentation"
        )
        return None

    # ------------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------------

    def infer(self, preprocessed: Any) -> Dict:
        """Run thyroid nodule detection and feature extraction.

        In mock mode, returns realistic nodule findings with TI-RADS
        features. In real mode, would run MONAI thyroid segmentation
        followed by nodule-level feature extraction.
        """
        if self.mock_mode:
            return self._mock_inference()

        logger.warning("Real inference not yet implemented; using mock results")
        return self._mock_inference()

    def _mock_inference(self) -> Dict:
        """Return realistic mock thyroid ultrasound findings.

        Simulates a thyroid ultrasound with two nodules:
            - A 22mm solid hypoechoic taller-than-wide nodule with punctate echogenic foci
              and irregular margins in the right lobe (TR5, 12 points)
            - A 14mm mixed cystic-solid isoechoic nodule with smooth
              margins in the left lobe (TR2, 2 points)

        AUDIT FIX (2026-11-03): this docstring said "TR5, 10 points" and "TR3, 2 points". The
        code computes 12 and TR2 (2 points maps to TR2 by this module's own documented chart).
        The figures are asserted in tests/test_radiology_reference_cases.py now, because this
        narrative is what reviewers read instead of running the code -- in the companion lung
        workflow the same drift is what hid a grading defect for a release.
        """
        return {
            "thyroid_volume_ml": 18.5,
            "isthmus_thickness_mm": 3.2,
            "background_parenchyma": "heterogeneous",
            "nodules": [
                {
                    "id": "nodule_1",
                    "location": "right lobe mid",
                    "laterality": "right",
                    "max_diameter_mm": 22.0,
                    "ap_diameter_mm": 18.5,
                    "transverse_diameter_mm": 20.1,
                    "longitudinal_diameter_mm": 22.0,
                    "composition": "solid",
                    "echogenicity": "hypoechoic",
                    "shape": "taller_than_wide",
                    "margin": "irregular",
                    "echogenic_foci": "punctate_echogenic_foci",
                    "vascularity": "increased_intranodular",
                    "detection_confidence": 0.96,
                },
                {
                    "id": "nodule_2",
                    "location": "left lobe inferior",
                    "laterality": "left",
                    "max_diameter_mm": 14.0,
                    "ap_diameter_mm": 10.2,
                    "transverse_diameter_mm": 12.8,
                    "longitudinal_diameter_mm": 14.0,
                    "composition": "mixed_cystic_solid",
                    "echogenicity": "isoechoic",
                    "shape": "wider_than_tall",
                    "margin": "smooth",
                    "echogenic_foci": "none",
                    "vascularity": "peripheral",
                    "detection_confidence": 0.91,
                },
            ],
        }

    # ------------------------------------------------------------------
    # Postprocessing
    # ------------------------------------------------------------------

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Calculate TI-RADS scores and determine FNA recommendations.

        Applies ACR TI-RADS (2017) point system:
            - Calculate points for each nodule from 5 feature categories
            - Map total points to TR level
            - Determine FNA recommendation based on TR level + size
            - Overall severity driven by highest-TR nodule

        Args:
            inference_result: Raw inference dict with 'nodules' list.

        Returns:
            WorkflowResult with TI-RADS classification and severity.
        """
        nodules = inference_result.get("nodules", [])

        if not nodules:
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.COMPLETED,
                findings=[{
                    "category": FindingCategory.NORMAL.value,
                    "description": "No thyroid nodules detected.",
                    "severity": FindingSeverity.NORMAL.value,
                    "tirads": TIRADS.TR1.value,
                }],
                measurements={
                    "nodule_count": 0.0,
                    "thyroid_volume_ml": inference_result.get("thyroid_volume_ml", 0.0),
                },
                classification=f"TI-RADS {TIRADS.TR1.value}",
                severity=FindingSeverity.NORMAL,
                nim_services_used=self.MODELS_USED,
            )

        findings = []
        overall_severity = FindingSeverity.NORMAL
        highest_tirads: Optional[TIRADS] = None
        ungraded: List[str] = []

        severity_order = [
            FindingSeverity.NORMAL,
            FindingSeverity.ROUTINE,
            FindingSeverity.SIGNIFICANT,
            FindingSeverity.URGENT,
            FindingSeverity.CRITICAL,
        ]

        for nodule in nodules:
            nodule_id = nodule.get("id", "unknown")
            location = nodule.get("location", "unspecified")
            laterality = nodule.get("laterality", "unspecified")
            # AUDIT FIX (2026-11-03): this defaulted to 0.0, which is below every published
            # FNA and follow-up threshold, so an unmeasured nodule was managed as the smallest
            # possible one. None now means "not measured" and is carried as such.
            raw_diameter = nodule.get("max_diameter_mm")
            max_diameter_mm = None if raw_diameter is None else float(raw_diameter)

            feature_score = score_tirads_features(nodule)

            features = {
                feature: nodule.get(feature) for feature in FEATURE_POINT_TABLES
            }
            measurements_block = {
                "max_diameter_mm": max_diameter_mm,
                "ap_diameter_mm": nodule.get("ap_diameter_mm"),
                "transverse_diameter_mm": nodule.get("transverse_diameter_mm"),
                "longitudinal_diameter_mm": nodule.get("longitudinal_diameter_mm"),
            }

            # --- A feature that could not be scored means no TR level at all -----------
            if not feature_score.is_complete:
                problems = feature_score.describe_problems()
                ungraded.append(f"{nodule_id} ({problems})")
                findings.append({
                    "category": FindingCategory.INDETERMINATE.value,
                    "description": (
                        f"Thyroid nodule in {location}: TI-RADS NOT ASSESSED -- {problems}. "
                        "An unscored ACR feature cannot be treated as a benign one, so no "
                        "TR level and no management recommendation are issued."
                    ),
                    "severity": FindingSeverity.UNKNOWN.value,
                    "nodule_id": nodule_id,
                    "laterality": laterality,
                    "location": location,
                    "tirads": None,
                    "tirads_points": None,
                    "unscored_features": {
                        "missing": sorted(feature_score.missing),
                        "unrecognised": dict(sorted(feature_score.unrecognised.items())),
                    },
                    "partial_points_scored": feature_score.points,
                    "recommendation": (
                        "TI-RADS cannot be assigned until every ACR feature is reported with "
                        "a label the published tables contain. Review the source study."
                    ),
                    "detection_confidence": nodule.get("detection_confidence", 0.0),
                    "features": features,
                    "measurements": measurements_block,
                })
                continue

            points = feature_score.points
            tr_level = points_to_tirads(points)
            nodule_severity = tirads_to_severity(tr_level)
            recommendation = tirads_fna_recommendation(tr_level, max_diameter_mm)

            if max_diameter_mm is None and tr_level in (
                TIRADS.TR3, TIRADS.TR4, TIRADS.TR5
            ):
                ungraded.append(f"{nodule_id} (max_diameter_mm not measured)")

            # Track overall severity (highest wins)
            if severity_order.index(nodule_severity) > severity_order.index(overall_severity):
                overall_severity = nodule_severity

            # AUDIT FIX (2026-11-03): highest_tirads is tracked by TR ordinal, not by the
            # severity index. TR1 and TR2 share one severity, so the old test could not see a
            # TR2 and the headline said TR1.
            if highest_tirads is None or TIRADS_ORDER[tr_level] > TIRADS_ORDER[highest_tirads]:
                highest_tirads = tr_level

            size_text = (
                f"{max_diameter_mm:.1f}mm" if max_diameter_mm is not None
                else "size not measured"
            )
            description = (
                f"Thyroid nodule in {location}, {size_text}, "
                f"{nodule.get('composition', 'unspecified')} composition, "
                f"{nodule.get('echogenicity', 'unspecified')} echogenicity, "
                f"TI-RADS {tr_level.value} ({points} points)"
            )

            findings.append({
                "category": FindingCategory.THYROID_NODULE.value,
                "description": description,
                "severity": nodule_severity.value,
                "nodule_id": nodule_id,
                "laterality": laterality,
                "location": location,
                "tirads": tr_level.value,
                "tirads_points": points,
                "recommendation": recommendation,
                "detection_confidence": nodule.get("detection_confidence", 0.0),
                "features": features,
                "measurements": measurements_block,
            })

        # Aggregate measurements. Only genuinely computed values go in here: a nodule whose
        # features could not be scored contributes no point total, because a 0.0 in this dict
        # is indistinguishable from a benign nodule.
        measurements: Dict[str, float] = {
            "nodule_count": float(len(nodules)),
            "thyroid_volume_ml": inference_result.get("thyroid_volume_ml", 0.0),
            "isthmus_thickness_mm": inference_result.get("isthmus_thickness_mm", 0.0),
        }
        for i, nodule in enumerate(nodules):
            prefix = f"nodule_{i + 1}"
            diameter = nodule.get("max_diameter_mm")
            if diameter is not None:
                measurements[f"{prefix}_diameter_mm"] = float(diameter)
            nodule_score = score_tirads_features(nodule)
            if nodule_score.is_complete:
                measurements[f"{prefix}_tirads_points"] = float(nodule_score.points)

        # Cross-modal trigger for TR4+ (BRAF, RAS, RET/PTC).
        #
        # AUDIT FIX (2026-11-03): this used to stamp findings[0] -- the FIRST nodule in the
        # list, not the nodule whose TR level fired the trigger. A gland with a benign simple
        # cyst listed first and a TR5 nodule second recorded the BRAF/RAS/RET-PTC work-up
        # against the cyst, in the wrong lobe, while the TR5 nodule carried no molecular
        # annotation. Anything joining genomics to a nodule ID -- or an FNA ordered from that
        # record -- then pointed at the wrong lesion. The flag now goes on each finding whose
        # OWN level is TR4 or TR5.
        triggering_findings = [
            f for f in findings
            if f.get("tirads") in (TIRADS.TR4.value, TIRADS.TR5.value)
        ]
        for finding in triggering_findings:
            finding["cross_modal_trigger"] = True
            finding["genomic_queries"] = ["BRAF V600E", "RAS", "RET/PTC"]

        # A gland in which nothing could be graded, or any nodule could not be graded, is not
        # a completed read. Say so in the status, the classification and the severity rather
        # than reporting the nodules that did grade and letting the rest look benign.
        if ungraded:
            reason = "; ".join(ungraded)
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.FAILED,
                findings=findings,
                measurements=measurements,
                classification=f"TI-RADS indeterminate -- {reason}",
                severity=FindingSeverity.UNKNOWN,
                nim_services_used=self.MODELS_USED,
            )

        assert highest_tirads is not None  # no ungraded nodules and nodules is non-empty

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=findings,
            measurements=measurements,
            classification=f"TI-RADS {highest_tirads.value}",
            severity=overall_severity,
            nim_services_used=self.MODELS_USED,
        )
