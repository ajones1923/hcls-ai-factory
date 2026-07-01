"""Evidence-based imaging protocol optimization.

Recommends optimal acquisition parameters based on clinical indication,
patient characteristics, and ACR Appropriateness Criteria.

Pure Python, no external licenses. Apache 2.0 compatible.

Author: Adam Jones
Date: April 2026
"""

from difflib import SequenceMatcher
from typing import Any, Dict, List, Optional

from loguru import logger
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════
# DATA MODELS
# ═══════════════════════════════════════════════════════════════════════


class PatientFactors(BaseModel):
    """Patient-specific factors that influence protocol selection."""

    age: Optional[int] = None
    weight_kg: Optional[float] = None
    bmi: Optional[float] = None
    sex: Optional[str] = None  # "M", "F"
    pregnant: bool = False
    renal_function_egfr: Optional[float] = None  # mL/min/1.73m2
    contrast_allergy: bool = False
    contrast_allergy_severity: Optional[str] = None  # "mild", "moderate", "severe"
    pediatric: bool = False


class ProtocolRecommendation(BaseModel):
    """Full protocol recommendation with evidence rating and patient adjustments."""

    indication: str = ""
    recommended_modality: str = ""
    recommended_protocol: str = ""
    acr_appropriateness_rating: int = Field(
        default=5, ge=1, le=9,
        description="1=usually not appropriate, 9=usually appropriate",
    )
    parameters: Dict[str, Any] = Field(default_factory=dict)
    contrast: Optional[Dict[str, Any]] = None
    dose_estimate_msv: Optional[float] = None
    warnings: List[str] = Field(default_factory=list)
    alternatives: List[Dict[str, Any]] = Field(default_factory=list)
    rationale: str = ""


# ═══════════════════════════════════════════════════════════════════════
# PROTOCOL OPTIMIZER
# ═══════════════════════════════════════════════════════════════════════


class ProtocolOptimizer:
    """Evidence-based protocol recommendation engine.

    Contains simplified ACR Appropriateness Criteria for 10+ common
    clinical indications, with patient-specific adjustments for
    pregnancy, renal impairment, contrast allergy, pediatric, and
    body habitus.
    """

    # ACR Appropriateness Criteria (simplified, key indications)
    ACR_CRITERIA: Dict[str, List[Dict[str, Any]]] = {
        "headache": [
            {
                "modality": "MRI",
                "protocol": "MRI Brain without contrast",
                "rating": 9,
                "contrast": None,
                "parameters": {"sequences": ["T1", "T2", "FLAIR", "DWI"], "slice_mm": 3.0},
            },
            {
                "modality": "CT",
                "protocol": "CT Head without contrast",
                "rating": 7,
                "contrast": None,
                "parameters": {"kv": 120, "mas": 250, "slice_mm": 5.0},
            },
            {
                "modality": "MRI",
                "protocol": "MRI Brain with contrast",
                "rating": 6,
                "contrast": {"agent": "gadolinium", "volume_ml": 15, "timing": "post-injection"},
                "parameters": {"sequences": ["T1", "T1+C", "FLAIR", "DWI"], "slice_mm": 3.0},
            },
        ],
        "chest_pain_acute": [
            {
                "modality": "CT",
                "protocol": "CTA Coronary",
                "rating": 9,
                "contrast": {"agent": "iodinated", "volume_ml": 80, "timing": "bolus-tracked"},
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 0.625, "ecg_gating": True},
            },
            {
                "modality": "CT",
                "protocol": "CT Chest PE protocol",
                "rating": 8,
                "contrast": {"agent": "iodinated", "volume_ml": 75, "timing": "bolus-tracked"},
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 1.25},
            },
            {
                "modality": "XR",
                "protocol": "Chest X-Ray PA/Lateral",
                "rating": 7,
                "contrast": None,
                "parameters": {"views": ["PA", "lateral"], "kv": 120},
            },
        ],
        "lung_cancer_screening": [
            {
                "modality": "CT",
                "protocol": "Low-dose CT Chest",
                "rating": 9,
                "contrast": None,
                "parameters": {"kv": 100, "mas": 30, "slice_mm": 1.0, "iterative_recon": True},
            },
            {
                "modality": "XR",
                "protocol": "Chest X-Ray",
                "rating": 2,
                "contrast": None,
                "parameters": {"views": ["PA"], "kv": 120},
            },
        ],
        "stroke_acute": [
            {
                "modality": "CT",
                "protocol": "CT Head without contrast",
                "rating": 9,
                "contrast": None,
                "parameters": {"kv": 120, "mas": 300, "slice_mm": 5.0},
            },
            {
                "modality": "CT",
                "protocol": "CTA Head/Neck",
                "rating": 9,
                "contrast": {"agent": "iodinated", "volume_ml": 80, "timing": "bolus-tracked"},
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 0.625},
            },
            {
                "modality": "MRI",
                "protocol": "MRI Brain stroke protocol",
                "rating": 8,
                "contrast": None,
                "parameters": {"sequences": ["DWI", "FLAIR", "MRA", "SWI"], "slice_mm": 3.0},
            },
        ],
        "abdominal_pain": [
            {
                "modality": "CT",
                "protocol": "CT Abdomen/Pelvis with contrast",
                "rating": 9,
                "contrast": {"agent": "iodinated", "volume_ml": 100, "timing": "portal venous phase"},
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 3.0},
            },
            {
                "modality": "US",
                "protocol": "Ultrasound Abdomen",
                "rating": 7,
                "contrast": None,
                "parameters": {"probe": "curvilinear", "frequency_mhz": "3-5"},
            },
            {
                "modality": "XR",
                "protocol": "Abdominal X-Ray",
                "rating": 5,
                "contrast": None,
                "parameters": {"views": ["supine", "upright"], "kv": 80},
            },
        ],
        "breast_screening": [
            {
                "modality": "MAMMO",
                "protocol": "Digital Mammography",
                "rating": 9,
                "contrast": None,
                "parameters": {"views": ["CC", "MLO"], "compression": True},
            },
            {
                "modality": "MRI",
                "protocol": "MRI Breast with contrast",
                "rating": 7,
                "contrast": {"agent": "gadolinium", "volume_ml": 15, "timing": "dynamic"},
                "parameters": {"sequences": ["T1", "T1+C", "STIR"], "slice_mm": 1.0},
                "note": "High-risk patients",
            },
            {
                "modality": "US",
                "protocol": "Breast Ultrasound",
                "rating": 5,
                "contrast": None,
                "parameters": {"probe": "linear", "frequency_mhz": "7-12"},
            },
        ],
        "pulmonary_embolism": [
            {
                "modality": "CT",
                "protocol": "CTA Chest PE protocol",
                "rating": 9,
                "contrast": {"agent": "iodinated", "volume_ml": 75, "timing": "bolus-tracked"},
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 1.25},
            },
            {
                "modality": "NM",
                "protocol": "V/Q Scan",
                "rating": 6,
                "contrast": {"agent": "Tc-99m MAA", "volume_ml": "IV injection", "timing": "perfusion"},
                "parameters": {"views": ["anterior", "posterior", "lateral", "oblique"]},
            },
        ],
        "liver_lesion": [
            {
                "modality": "MRI",
                "protocol": "MRI Abdomen liver protocol",
                "rating": 9,
                "contrast": {"agent": "gadolinium (hepatobiliary)", "volume_ml": 15, "timing": "multiphase"},
                "parameters": {"sequences": ["T1", "T2", "DWI", "dynamic"], "slice_mm": 3.0},
            },
            {
                "modality": "CT",
                "protocol": "CT Abdomen multiphase",
                "rating": 8,
                "contrast": {"agent": "iodinated", "volume_ml": 120, "timing": "arterial + portal venous + delayed"},
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 3.0},
            },
            {
                "modality": "US",
                "protocol": "Ultrasound Abdomen",
                "rating": 6,
                "contrast": None,
                "parameters": {"probe": "curvilinear", "frequency_mhz": "3-5"},
            },
        ],
        "thyroid_nodule": [
            {
                "modality": "US",
                "protocol": "Thyroid Ultrasound",
                "rating": 9,
                "contrast": None,
                "parameters": {"probe": "linear", "frequency_mhz": "7-15"},
            },
            {
                "modality": "NM",
                "protocol": "Thyroid Scintigraphy",
                "rating": 5,
                "contrast": {"agent": "I-123 or Tc-99m", "volume_ml": "oral/IV", "timing": "4-24 hour uptake"},
                "parameters": {"views": ["anterior"]},
            },
        ],
        "prostate_cancer": [
            {
                "modality": "MRI",
                "protocol": "mpMRI Prostate (PI-RADS)",
                "rating": 9,
                "contrast": {"agent": "gadolinium", "volume_ml": 15, "timing": "dynamic"},
                "parameters": {"sequences": ["T2", "DWI", "DCE"], "slice_mm": 3.0},
            },
            {
                "modality": "US",
                "protocol": "TRUS Prostate",
                "rating": 6,
                "contrast": None,
                "parameters": {"probe": "endorectal", "frequency_mhz": "7-9"},
            },
        ],
        "kidney_stone": [
            {
                "modality": "CT",
                "protocol": "CT Abdomen/Pelvis without contrast",
                "rating": 9,
                "contrast": None,
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 3.0, "low_dose": True},
            },
            {
                "modality": "US",
                "protocol": "Renal Ultrasound",
                "rating": 7,
                "contrast": None,
                "parameters": {"probe": "curvilinear", "frequency_mhz": "3-5"},
            },
        ],
        "appendicitis": [
            {
                "modality": "CT",
                "protocol": "CT Abdomen/Pelvis with contrast",
                "rating": 9,
                "contrast": {"agent": "iodinated", "volume_ml": 100, "timing": "portal venous phase"},
                "parameters": {"kv": 120, "mas": "auto", "slice_mm": 3.0},
            },
            {
                "modality": "US",
                "protocol": "Ultrasound Abdomen RLQ",
                "rating": 8,
                "contrast": None,
                "parameters": {"probe": "linear + curvilinear", "frequency_mhz": "5-12"},
                "note": "Preferred first-line in pediatric and pregnant patients",
            },
        ],
    }

    # CT protocol dose estimates (mSv effective dose)
    DOSE_ESTIMATES: Dict[str, float] = {
        "CT Head without contrast": 2.0,
        "CTA Head/Neck": 5.0,
        "CT Chest PE protocol": 7.0,
        "CTA Chest PE protocol": 7.0,
        "Low-dose CT Chest": 1.5,
        "CT Abdomen/Pelvis with contrast": 10.0,
        "CT Abdomen/Pelvis without contrast": 8.0,
        "CT Abdomen multiphase": 15.0,
        "CTA Coronary": 5.0,
        "Chest X-Ray PA/Lateral": 0.1,
        "Chest X-Ray": 0.02,
        "Abdominal X-Ray": 0.7,
        "Digital Mammography": 0.4,
        "V/Q Scan": 2.0,
        "Thyroid Scintigraphy": 1.0,
    }

    # Contrast agents by type
    IODINATED_CONTRAST_MODALITIES = {"CT"}
    GADOLINIUM_CONTRAST_MODALITIES = {"MRI"}

    # Keyword synonyms for fuzzy matching
    _INDICATION_SYNONYMS: Dict[str, List[str]] = {
        "headache": ["headache", "migraine", "cephalgia", "head pain"],
        "chest_pain_acute": ["chest pain", "angina", "acs", "cardiac pain", "mi"],
        "lung_cancer_screening": ["lung screening", "lung cancer", "ldct", "low dose ct lung"],
        "stroke_acute": ["stroke", "cva", "tia", "cerebrovascular", "ischemic stroke", "brain attack"],
        "abdominal_pain": ["abdominal pain", "belly pain", "stomach pain", "acute abdomen"],
        "breast_screening": ["breast screening", "mammogram", "breast cancer screening"],
        "pulmonary_embolism": ["pulmonary embolism", "pe", "dvt", "venous thromboembolism"],
        "liver_lesion": ["liver lesion", "hepatic mass", "liver mass", "hcc", "liver cancer"],
        "thyroid_nodule": ["thyroid nodule", "thyroid mass", "goiter"],
        "prostate_cancer": ["prostate cancer", "prostate", "psa elevated", "prostate mass"],
        "kidney_stone": ["kidney stone", "renal calculus", "urolithiasis", "nephrolithiasis", "renal colic"],
        "appendicitis": ["appendicitis", "rlq pain", "right lower quadrant"],
    }

    def __init__(self) -> None:
        logger.info("ProtocolOptimizer initialized with {} indications", len(self.ACR_CRITERIA))

    # ───────────────────────────────────────────────────────────────────
    # PUBLIC API
    # ───────────────────────────────────────────────────────────────────

    def recommend(
        self,
        indication: str,
        patient: Optional[PatientFactors] = None,
    ) -> ProtocolRecommendation:
        """Recommend optimal imaging protocol based on indication and patient factors.

        Parameters
        ----------
        indication : str
            Clinical indication (e.g. "headache", "chest pain", "lung screening").
        patient : PatientFactors, optional
            Patient-specific factors for protocol adjustment.

        Returns
        -------
        ProtocolRecommendation
            Full recommendation including modality, protocol, ACR rating,
            dose estimate, warnings, and alternatives.
        """
        matched_key = self._match_indication(indication)
        logger.info("Indication '{}' matched to '{}'", indication, matched_key)

        criteria = self.ACR_CRITERIA.get(matched_key, [])
        if not criteria:
            logger.warning("No ACR criteria found for '{}', returning default", matched_key)
            return ProtocolRecommendation(
                indication=indication,
                recommended_modality="CT",
                recommended_protocol="CT per clinical judgment",
                acr_appropriateness_rating=5,
                parameters={},
                rationale=f"No specific ACR criteria matched for '{indication}'. "
                          "Clinical judgment recommended.",
            )

        best = criteria[0]

        rec = ProtocolRecommendation(
            indication=indication,
            recommended_modality=best["modality"],
            recommended_protocol=best["protocol"],
            acr_appropriateness_rating=best["rating"],
            parameters=dict(best.get("parameters", {})),
            contrast=dict(best["contrast"]) if best.get("contrast") else None,
            dose_estimate_msv=self.DOSE_ESTIMATES.get(best["protocol"]),
            alternatives=[
                {
                    "modality": alt["modality"],
                    "protocol": alt["protocol"],
                    "rating": alt["rating"],
                    "note": alt.get("note", ""),
                }
                for alt in criteria[1:]
            ],
            rationale=self._build_rationale(matched_key, best),
        )

        # Apply patient-specific adjustments
        if patient:
            rec = self._adjust_for_patient(rec, patient)

        logger.info(
            "Recommendation: {} — {} (rating {})",
            rec.recommended_modality,
            rec.recommended_protocol,
            rec.acr_appropriateness_rating,
        )
        return rec

    def get_available_indications(self) -> List[str]:
        """Return all supported clinical indications."""
        return sorted(self.ACR_CRITERIA.keys())

    # ───────────────────────────────────────────────────────────────────
    # INDICATION MATCHING
    # ───────────────────────────────────────────────────────────────────

    def _match_indication(self, indication: str) -> str:
        """Fuzzy-match user indication to ACR_CRITERIA keys.

        Uses synonym lists and SequenceMatcher for best match.
        """
        indication_lower = indication.lower().strip()

        # Exact key match
        if indication_lower in self.ACR_CRITERIA:
            return indication_lower

        # Synonym matching — check if any synonym appears in the indication
        best_key = ""
        best_score = 0.0

        for key, synonyms in self._INDICATION_SYNONYMS.items():
            for synonym in synonyms:
                if synonym in indication_lower:
                    score = len(synonym) / max(len(indication_lower), 1)
                    if score > best_score:
                        best_score = score
                        best_key = key

        if best_key and best_score > 0.15:
            return best_key

        # Fallback: SequenceMatcher against all keys
        for key in self.ACR_CRITERIA:
            ratio = SequenceMatcher(None, indication_lower, key).ratio()
            if ratio > best_score:
                best_score = ratio
                best_key = key

        if best_key and best_score > 0.4:
            return best_key

        # No match — return the raw indication for default handling
        logger.warning("No indication match for '{}' (best={}, score={:.2f})", indication, best_key, best_score)
        return ""

    # ───────────────────────────────────────────────────────────────────
    # PATIENT ADJUSTMENTS
    # ───────────────────────────────────────────────────────────────────

    def _adjust_for_patient(
        self,
        rec: ProtocolRecommendation,
        patient: PatientFactors,
    ) -> ProtocolRecommendation:
        """Apply patient-specific adjustments to the recommendation."""
        warnings: List[str] = list(rec.warnings)

        # Pregnancy adjustments
        if patient.pregnant:
            warnings.extend(self._adjust_for_pregnancy(rec))

        # Renal function adjustments
        if patient.renal_function_egfr is not None:
            warnings.extend(self._check_contrast_safety(patient))

        # Contrast allergy adjustments
        if patient.contrast_allergy:
            warnings.extend(self._adjust_for_contrast_allergy(rec, patient))

        # Commit accumulated warnings before pediatric (which returns new rec)
        rec = rec.model_copy(update={"warnings": warnings})

        # Pediatric adjustments (modifies rec including warnings, params, dose)
        if patient.pediatric or (patient.age is not None and patient.age < 18):
            rec = self._adjust_pediatric(rec, patient)

        # BMI / body habitus adjustments
        if patient.bmi is not None and patient.bmi > 35:
            updated_warnings = list(rec.warnings)
            updated_warnings.append(
                "Large body habitus (BMI > 35): consider increasing kV to 140, "
                "increase mAs, use larger FOV. Image noise may be elevated."
            )
            params = dict(rec.parameters)
            if "kv" in params and isinstance(params["kv"], int) and params["kv"] < 140:
                params["kv"] = 140
            rec = rec.model_copy(update={"parameters": params, "warnings": updated_warnings})
        elif patient.weight_kg is not None and patient.weight_kg > 120:
            updated_warnings = list(rec.warnings)
            updated_warnings.append(
                "Patient weight > 120 kg: consider increasing kV/mAs for adequate penetration."
            )
            rec = rec.model_copy(update={"warnings": updated_warnings})

        return rec

    def _adjust_for_pregnancy(self, rec: ProtocolRecommendation) -> List[str]:
        """Generate warnings for pregnant patients."""
        warnings: List[str] = []
        ionizing = {"CT", "XR", "NM", "MAMMO", "FLUORO"}

        if rec.recommended_modality.upper() in ionizing:
            warnings.append(
                f"PREGNANCY WARNING: {rec.recommended_modality} involves ionizing radiation. "
                "Prefer MRI or ultrasound if clinically appropriate. "
                "If CT is essential, minimize scan range and use lowest possible dose."
            )
        if rec.contrast and rec.contrast.get("agent", "").lower().startswith("iod"):
            warnings.append(
                "PREGNANCY WARNING: Iodinated contrast crosses the placenta. "
                "Use only if benefit outweighs risk. Neonatal thyroid screening recommended."
            )
        if rec.contrast and "gadolinium" in rec.contrast.get("agent", "").lower():
            warnings.append(
                "PREGNANCY WARNING: Gadolinium is classified as pregnancy risk category C. "
                "Use only if essential. Macrocyclic agents preferred."
            )
        return warnings

    def _check_contrast_safety(self, patient: PatientFactors) -> List[str]:
        """Check renal function for contrast safety."""
        warnings: List[str] = []
        egfr = patient.renal_function_egfr

        if egfr is not None and egfr < 30:
            warnings.append(
                f"RENAL IMPAIRMENT (eGFR {egfr} mL/min): "
                "Iodinated contrast is relatively contraindicated. "
                "Risk of contrast-induced nephropathy. "
                "Consider non-contrast study or alternative modality."
            )
            warnings.append(
                f"RENAL IMPAIRMENT (eGFR {egfr} mL/min): "
                "Gadolinium-based contrast agents carry risk of nephrogenic systemic "
                "fibrosis (NSF). Use Group II agents only if essential."
            )
        elif egfr is not None and egfr < 60:
            warnings.append(
                f"RENAL CAUTION (eGFR {egfr} mL/min): "
                "Moderate renal impairment. Pre/post hydration recommended if "
                "iodinated contrast is used. Monitor creatinine 48-72h post-procedure."
            )

        return warnings

    def _adjust_for_contrast_allergy(
        self,
        rec: ProtocolRecommendation,
        patient: PatientFactors,
    ) -> List[str]:
        """Generate warnings and adjustments for contrast allergy."""
        warnings: List[str] = []
        severity = (patient.contrast_allergy_severity or "unknown").lower()

        if rec.contrast:
            if severity == "severe":
                warnings.append(
                    "SEVERE CONTRAST ALLERGY: Contrast is contraindicated. "
                    "Use non-contrast protocol or alternative modality. "
                    "If contrast is absolutely necessary, full premedication "
                    "protocol required with anesthesia standby."
                )
            elif severity == "moderate":
                warnings.append(
                    "MODERATE CONTRAST ALLERGY: Premedication required "
                    "(prednisone 50 mg PO at 13h, 7h, and 1h prior; "
                    "diphenhydramine 50 mg IV/PO 1h prior). "
                    "Use non-ionic low-osmolarity contrast."
                )
            else:
                warnings.append(
                    "CONTRAST ALLERGY REPORTED: Premedication recommended. "
                    "Use non-ionic low-osmolarity contrast agent."
                )
        return warnings

    def _adjust_pediatric(
        self,
        rec: ProtocolRecommendation,
        patient: PatientFactors,
    ) -> ProtocolRecommendation:
        """Adjust recommendation for pediatric patients."""
        warnings = list(rec.warnings)
        params = dict(rec.parameters)

        warnings.append(
            "PEDIATRIC PATIENT: Use age/weight-appropriate protocols. "
            "Apply ALARA principle. Use size-specific dose estimates (SSDE)."
        )

        # Reduce kV for CT
        if rec.recommended_modality.upper() == "CT":
            if "kv" in params and isinstance(params["kv"], int):
                original_kv = params["kv"]
                if patient.age is not None and patient.age < 5:
                    params["kv"] = min(params["kv"], 80)
                elif patient.age is not None and patient.age < 12:
                    params["kv"] = min(params["kv"], 100)
                else:
                    params["kv"] = min(params["kv"], 120)

                if params["kv"] < original_kv:
                    warnings.append(
                        f"Pediatric dose reduction: kV reduced from {original_kv} to {params['kv']}."
                    )

            # Append pediatric protocol note
            protocol = rec.recommended_protocol
            if "pediatric" not in protocol.lower():
                protocol = f"{protocol} (Pediatric)"

            # Reduce dose estimate
            dose = rec.dose_estimate_msv
            if dose is not None:
                age = patient.age or 10
                if age < 1:
                    dose *= 0.3
                elif age < 5:
                    dose *= 0.5
                elif age < 10:
                    dose *= 0.6
                elif age < 15:
                    dose *= 0.8

            return rec.model_copy(update={
                "warnings": warnings,
                "parameters": params,
                "recommended_protocol": protocol,
                "dose_estimate_msv": round(dose, 2) if dose is not None else None,
            })

        return rec.model_copy(update={"warnings": warnings, "parameters": params})

    # ───────────────────────────────────────────────────────────────────
    # HELPERS
    # ───────────────────────────────────────────────────────────────────

    def _build_rationale(self, matched_key: str, best: Dict[str, Any]) -> str:
        """Build human-readable rationale for the recommendation."""
        rating = best["rating"]
        modality = best["modality"]
        protocol = best["protocol"]

        if rating >= 8:
            appropriateness = "usually appropriate"
        elif rating >= 5:
            appropriateness = "may be appropriate"
        else:
            appropriateness = "usually not appropriate"

        rationale = (
            f"Based on ACR Appropriateness Criteria for '{matched_key}', "
            f"{protocol} ({modality}) is rated {rating}/9 ({appropriateness}). "
        )

        if best.get("note"):
            rationale += f"Note: {best['note']}. "

        if best.get("contrast"):
            agent = best["contrast"].get("agent", "")
            rationale += f"Contrast agent: {agent}. "

        return rationale.strip()
